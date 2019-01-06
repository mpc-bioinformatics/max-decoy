extern crate postgres;
extern crate postgres_array;

use std::hash::{Hash, Hasher};
use std::collections::HashMap;
use std::collections::hash_map::Entry::{Occupied, Vacant};

use self::postgres::Connection;
use self::postgres_array::Array;

use proteomic::models::decoys::decoy::{Decoy, PlainDecoy};
use proteomic::models::decoys::base_decoy::BaseDecoy;
use proteomic::models::persistable::{handle_postgres_error, Persistable, QueryError, QueryOk, FromSqlRowError};
use proteomic::models::amino_acids::modification::Modification;
use proteomic::models::mass;
use proteomic::utility::database_connection::DatabaseConnection;

pub struct ModifiedDecoy {
    id: i64,                                        // SERIAL
    base_decoy_id: i64,                             // BIGINT REFERENCE
    base_decoy: BaseDecoy,                          
    c_terminus_modification: Option<Modification>,  
    c_terminus_modification_id: i64,                // BIGINT REFERENCE
    n_terminus_modification: Option<Modification>,
    n_terminus_modification_id: i64,                // BIGINT REFERENCE
    modifications: Vec<Option<Modification>>,
    modification_ids: Array<i64>,                   // BIGINT ARRAY
    weight: i64,                                    // BIGINT
    header_addition: String                         // TEXT
}

impl ModifiedDecoy {
    pub fn new(base_decoy: &BaseDecoy, c_terminus_modification: Option<Modification>, n_terminus_modification: Option<Modification>, modifications: &Vec<Option<Modification>>, weight: i64) -> ModifiedDecoy {
        let dummy_modification = Modification::get_dummy_modification();
        let header_addition = Self::create_header_addition(&c_terminus_modification, &n_terminus_modification, modifications);
        return Self {
            id: 0,
            base_decoy_id: base_decoy.get_primary_key(),
            base_decoy: base_decoy.clone(),
            c_terminus_modification_id: match &c_terminus_modification {
                Some(modification) => modification.get_primary_key(),
                None => dummy_modification.get_primary_key()
            },
            c_terminus_modification: c_terminus_modification,
            n_terminus_modification_id: match &n_terminus_modification { 
                Some(modification) => modification.get_primary_key(),
                None => dummy_modification.get_primary_key()
            },
            n_terminus_modification: n_terminus_modification,
            modification_ids: *Self::modifications_to_psql_array(modifications),
            modifications: modifications.clone(),
            weight: weight,
            header_addition: header_addition
        }
    }

    fn create_header_addition(c_terminus_modification: &Option<Modification>, n_terminus_modification: &Option<Modification>, modifications: &Vec<Option<Modification>>) -> String {
    // create a HashMap with key "mod_accession|mod_name" and value Vec<String> which contains the positions of the modification
        let mut accession_header_position_map: HashMap<String, Vec<String>> = HashMap::new();
        for idx in 0..(modifications.len()) {
            match &modifications[idx] {
                Some(modification) => {
                    let peff_notification_of_accession_and_name = format!("{}|{}", modification.get_accession(), modification.get_name());
                    let positions = match accession_header_position_map.entry(peff_notification_of_accession_and_name) {
                        Vacant(entry) => entry.insert(Vec::new()),
                        Occupied(entry) => entry.into_mut(),
                    };
                    positions.push(idx.to_string());
                },
                None => ()
            }
        }
        match c_terminus_modification {
            Some(modification) => {
                let peff_notification_of_accession_and_name = format!("{}|{}", modification.get_accession(), modification.get_name());
                let positions = match accession_header_position_map.entry(peff_notification_of_accession_and_name) {
                    Vacant(entry) => entry.insert(Vec::new()),
                    Occupied(entry) => entry.into_mut(),
                };
                positions.push("c-terminus".to_string());
            },
            None => ()
        }
        match n_terminus_modification {
            Some(modification) => {
                let peff_notification_of_accession_and_name = format!("{}|{}", modification.get_accession(), modification.get_name());
                let positions = match accession_header_position_map.entry(peff_notification_of_accession_and_name) {
                    Vacant(entry) => entry.insert(Vec::new()),
                    Occupied(entry) => entry.into_mut(),
                };
                positions.push("n-terminus".to_string());
            },
            None => ()
        }
        let mut mod_res: String = String::new();
        for (accession_and_header, positions) in &accession_header_position_map {
            mod_res.push_str(
                format!(
                    "({}|{})",
                    positions.join(","),
                    accession_and_header
                ).as_str()
            );
        }
        return format!("ModRes={}", mod_res);
    }

    fn modifications_to_psql_array(modifications: &Vec<Option<Modification>>) -> Box<Array<i64>> {
        let dummy_modification = Modification::get_dummy_modification();
        let mut modification_ids: Vec<i64> = Vec::new();
        for modification_option in modifications.iter() {
            match modification_option {
                Some(modification) => modification_ids.push(modification.get_primary_key()),
                None => modification_ids.push(dummy_modification.get_primary_key())
            }
        }
        return Box::new(Array::from_vec(modification_ids, modifications.len() as i32));
    }

    fn modifications_from_psql_array(conn: &postgres::Connection, modification_ids: &Array<i64>) -> Result<Vec<Option<Modification>>, QueryError> {
        let dummy_modification = Modification::get_dummy_modification();
        // cast ids to string and join them with ', '
        let modification_ids_for_query: String = modification_ids.iter().map(|id| id.to_string()).collect::<Vec<String>>().join(", ");
        // condition for 'SELECT amino_acid_modifications WHERE id IN (id1, id2, ...)' 
        let condition: String = format!("id IN ({})", modification_ids_for_query);
        let modifications = match Modification::find_where(conn, condition.as_str(), &[]) {
            Ok(modifications) => modifications,
            Err(err) => return Err(err)
        };
        // check if for every id a modification is found
        if modifications.len() != modification_ids.iter().count() {
            // collect ids of missing modifications
            let mut ids_of_missing_modifications: Vec<i64> = Vec::new();
            'id_loop: for id in modification_ids.iter() {
                'mod_loop: for modification in modifications.iter() {
                    if modification.get_primary_key() == *id {
                        continue 'id_loop;
                    }
                }
                ids_of_missing_modifications.push(*id);
            }
            let missing_ids_string: String = ids_of_missing_modifications.iter().map(|id| id.to_string()).collect::<Vec<String>>().join(", ");
            return Err(QueryError::AssociatedRecordNotFound(format!("Modifications(ids: [{}])", missing_ids_string)));
        }
        let mut modification_options: Vec<Option<Modification>> = Vec::new();
        // wrap modifications in Option
        for modification in modifications.iter(){
            // replace all dummy-modifications with None
            if modification.get_primary_key() == dummy_modification.get_primary_key(){
                modification_options.push(None);
            } else {
                modification_options.push(Some(modification.clone()));
            }
        }
        return Ok(modification_options);
    }

    pub fn get_c_terminus_modification(&self) -> &Option<Modification> {
        return &self.c_terminus_modification;
    }

    pub fn get_n_terminus_modification(&self) -> &Option<Modification> {
        return &self.n_terminus_modification
    }

    fn get_c_terminus_modification_id(&self) -> &i64 {
        return &self.c_terminus_modification_id;
    }

    fn get_n_terminus_modification_id(&self) -> &i64 {
        return &self.n_terminus_modification_id;
    }

    fn get_modification_ids(&self) -> &Array<i64> {
        return &self.modification_ids;
    }

    // row must contain [base_decoy.aa_sequence, base_decoy.header, modified_decoy.header_addition, modified_decoy.weight] in this order
    fn prepared_decoy_from_sql_row(row: &postgres::rows::Row) -> PlainDecoy {
        return PlainDecoy::new(
            format!("{} {}", row.get::<usize, String>(1), row.get::<usize, String>(2)).as_str(),
            &row.get::<usize, String>(0), 
            &row.get(3)
        );
    }

    // this function works like find_where() but to save time it will return PlainDecoy instead of loading all modifications first
    pub fn find_where_as_prepared_decoys(conn: &postgres::Connection, conditions: &str, values: &[&postgres::types::ToSql]) -> Result<Vec<PlainDecoy>, QueryError> {
        let select_query: String = format!(
            "SELECT {base_decoys_table}.aa_sequence, {base_decoys_table}.header, {modified_decoys_table}.header_addition, {modified_decoys_table}.weight FROM {modified_decoys_table} INNER JOIN base_decoys ON {modified_decoys_table}.base_decoy_id={base_decoys_table}.id WHERE {conditions};",
            modified_decoys_table = Self::get_table_name(),
            base_decoys_table = BaseDecoy::get_table_name(),
            conditions = conditions
        );
        println!("{}", select_query);
        match conn.query(select_query.as_str(), values) {
            Ok(ref rows) => {
                let mut records: Vec<PlainDecoy> = Vec::new();
                for row in rows {
                    records.push(Self::prepared_decoy_from_sql_row(&row));
                }
                return Ok(records);
            },
            Err(err) => Err(handle_postgres_error(&err))
        }
    }
}

impl Decoy for ModifiedDecoy {
    fn to_string(&self) -> String {
        return format!(
            "proteomic::modes::decoys::modified_decoy::ModifiedDecoy\n\theader => {}\n\taa_sequence => {}\n\tweight => {}",
            self.get_header(),
            self.get_aa_sequence(),
            mass::convert_mass_to_float(self.get_weight())
        );
    }

    fn get_header(&self) -> String {
        return format!(
            "{} {}",
            self.base_decoy.get_header(),
            self.header_addition
        )
    }

    fn get_aa_sequence(&self) -> String {
        return self.base_decoy.get_aa_sequence();
    }

    fn get_weight(&self) -> i64 {
        return self.weight;
    }

    fn get_length(&self) -> i32 {
        return self.base_decoy.get_length() as i32;
    }

    fn get_c_terminus_amino_acid(&self) -> char {
        return self.base_decoy.get_c_terminus_amino_acid();
    }

    fn get_n_terminus_amino_acid(&self) -> char {
        return self.base_decoy.get_n_terminus_amino_acid();
    }

    fn get_amino_acid_at(&self, idx: usize) -> char {
        return self.base_decoy.get_amino_acid_at(idx);
    }
}


impl Persistable<ModifiedDecoy, i64, (i64, i64, i64, &Array<i64>, i64)> for ModifiedDecoy {
    fn from_sql_row(row: &postgres::rows::Row) -> Result<Self, FromSqlRowError> {
        let conn = DatabaseConnection::get_database_connection();
        let base_decoy: BaseDecoy = match BaseDecoy::find(&conn, &[&row.get::<usize, i64>(1)]) {
            Ok(base_decoy) => base_decoy,
            Err(query_error) => match query_error {
                QueryError::NoMatch => return Err(FromSqlRowError::AssociatedRecordNotFound(format!("BaseDecoy(id: {})", row.get::<usize, i64>(1)))),
                _ => return Err(FromSqlRowError::InnerQueryError(query_error))                          // else a sql error occures
            }
        };

        let dummy_modification: Modification = Modification::get_dummy_modification();

        let c_terminus_modification_id: i64 = row.get::<usize, i64>(2);
        let n_terminus_modification_id: i64 = row.get::<usize, i64>(3);

        let mut c_terminus_modification: Option<Modification> = None;
        if c_terminus_modification_id != dummy_modification.get_primary_key() {
            match Modification::find(&conn, &[&c_terminus_modification_id]) {
                Ok(modification) => c_terminus_modification = Some(modification),
                Err(error) => match error {
                    QueryError::NoMatch => return Err(FromSqlRowError::AssociatedRecordNotFound(format!("CTerminusModification(id: {})", c_terminus_modification_id))),
                    _ => return Err(FromSqlRowError::InnerQueryError(error))  // else a sql error occures
                }
            }
        }
            
        let mut n_terminus_modification: Option<Modification> = None;
        if n_terminus_modification_id != dummy_modification.get_primary_key() {
            match Modification::find(&conn, &[&n_terminus_modification_id]) {
                Ok(modification) => n_terminus_modification = Some(modification),
                Err(error) => match error {
                    QueryError::NoMatch => return Err(FromSqlRowError::AssociatedRecordNotFound(format!("NTerminusModification(id: {})", n_terminus_modification_id))),
                    _ => return Err(FromSqlRowError::InnerQueryError(error))  // else a sql error occures
                }
            }
        }

        let modification_ids: Array<i64> = row.get::<usize, Array<i64>>(4);
        let modifications: Vec<Option<Modification>> = match Self::modifications_from_psql_array(&conn, &modification_ids) {
            Ok(modifications) => modifications,
            Err(err) => return Err(FromSqlRowError::InnerQueryError(err))
        };

        return Ok (
            Self {
                id: row.get(0),
                base_decoy_id: base_decoy.get_primary_key(),
                base_decoy: base_decoy,
                c_terminus_modification_id: c_terminus_modification_id,
                c_terminus_modification: c_terminus_modification,
                n_terminus_modification_id: n_terminus_modification_id,
                n_terminus_modification: n_terminus_modification,
                modification_ids: modification_ids,
                modifications: modifications,
                weight: row.get(5),
                header_addition: row.get(6)
            }
        );
    }

    fn set_primary_key_from_sql_row(&mut self, row: &postgres::rows::Row) {
        self.id = row.get(0);
    }

    fn invalidate_primary_key(&mut self) {
        self.id = 0;
    }

    fn get_primary_key(&self) -> i64 {
        return self.id;
    }

    fn get_table_name() -> &'static str {
        return "modified_decoys";
    }

    fn is_persisted(&self) -> bool {
        return self.id > 0;
    }


    fn find_query() -> &'static str {
        return "SELECT * FROM modified_decoys WHERE id = $1 LIMIT 1;";
    }


    fn create_query() -> &'static str {
        return "INSERT INTO modified_decoys (base_decoy_id, c_terminus_modification_id, n_terminus_modification_id, modification_ids, weight, header_addition) VALUES ($1, $2, $3, $4, $5, $6) ON CONFLICT (base_decoy_id, c_terminus_modification_id, n_terminus_modification_id, modification_ids, weight) DO NOTHING RETURNING id;";
    }

    fn create_attributes(&self) -> Box<Vec<&postgres::types::ToSql>>{
        return Box::new(vec![&self.base_decoy_id, &self.c_terminus_modification_id, &self.n_terminus_modification_id, &self.modification_ids, &self.weight, &self.header_addition]);
    }

    fn update_query() -> &'static str{
        return "UPDATE modified_decoys SET base_decoy_id = $2, c_terminus_modification_id = $3, n_terminus_modification_id = $4, modification_ids = $5, weight = $6, header_addition = $7 WHERE id = $1;";
    }

    fn update_attributes(&self) -> Box<Vec<&postgres::types::ToSql>>{
        return Box::new(vec![&self.id, &self.base_decoy_id, &self.c_terminus_modification_id, &self.n_terminus_modification_id, &self.modification_ids, &self.weight, &self.header_addition]);
    }

    fn delete_query() -> &'static str {
        return "DELETE FROM modified_decoys WHERE id = $1;";
    }

    fn delete_attributes(&self) -> Box<Vec<&postgres::types::ToSql>> {
        return Box::new(vec![&self.id]);
    }

    fn delete_all_query() -> &'static str {
        return "DELETE FROM modified_decoys WHERE id IS NOT NULL;";
    }

    fn exists_query() -> &'static str {
        return "SELECT id FROM modified_decoys WHERE base_decoy_id = $1 AND c_terminus_modification_id = $2 AND n_terminus_modification_id = $3 AND modification_ids = $4 AND weight = $5;"
    }

    fn exists_attributes(&self) -> Box<Vec<&postgres::types::ToSql>> {
        return Box::new(vec![&self.base_decoy_id, &self.c_terminus_modification_id, &self.n_terminus_modification_id, &self.modification_ids, &self.weight])
    }

    fn before_delete_hook(&self) -> Result<(), QueryError> {return Ok(());}
}

// PartialEq-implementation to use this type in a HashSet
impl PartialEq for ModifiedDecoy {
    fn eq(&self, other: &ModifiedDecoy) -> bool {
       return (self.get_aa_sequence() == *other.get_aa_sequence()) 
        & (self.weight == other.get_weight())
        & (self.get_header() == *other.get_header());
    }
}

// Eq-implementation to use this type in a HashSet
impl Eq for ModifiedDecoy {}

// Hash-implementation to use this type in a HashSet
impl Hash for ModifiedDecoy {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.get_aa_sequence().hash(state);
        self.get_header().hash(state);
    }
}

impl Clone for ModifiedDecoy {
    fn clone(&self) -> Self {
        return Self{
            id: self.id,
            base_decoy_id: self.base_decoy_id,
            base_decoy: self.base_decoy.clone(),
            c_terminus_modification_id: self.c_terminus_modification_id,
            c_terminus_modification: self.c_terminus_modification.clone(),
            n_terminus_modification_id: self.n_terminus_modification_id,
            n_terminus_modification: self.n_terminus_modification.clone(),
            modification_ids: self.modification_ids.clone(),
            modifications: self.modifications.clone(),
            weight: self.get_weight(),
            header_addition: self.header_addition.clone()
        }
    }
}