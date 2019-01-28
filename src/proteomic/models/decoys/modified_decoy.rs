use std::hash::{Hash, Hasher};
use std::collections::HashMap;
use std::collections::hash_map::Entry::{Occupied, Vacant};

use postgres_array::Array;

use proteomic::models::decoys::decoy::Decoy;
use proteomic::models::decoys::base_decoy::BaseDecoy;
use proteomic::models::persistable::{handle_postgres_error, Persistable, QueryError, FromSqlRowError};
use proteomic::models::amino_acids::modification::Modification;
use proteomic::models::mass;
use proteomic::models::peptides::peptide_interface::PeptideInterface;
use proteomic::utility::database_connection::DatabaseConnection;

pub struct ModifiedDecoy {
    id: i64,                                        // SERIAL, primary key
    base_decoy_id: i64,                             // BIGINT REFERENCE, part of UNIQUE-constraints
    base_decoy: BaseDecoy,
    modifications: Vec<Modification>,
    modification_ids: Array<i64>,                   // BIGINT ARRAY, indexed
    weight: i64,                                    // BIGINT
    modification_summary: String,                   // TEXT, part of UNIQUE-constraints
    header: String
}

impl ModifiedDecoy {
    pub fn new(base_decoy: &BaseDecoy, modifications: &Vec<Modification>, weight: i64) -> ModifiedDecoy {
        let modification_summary = Self::create_modification_summary(modifications);
        return Self {
            id: 0,
            header: format!("{} ModSum={}", base_decoy.get_header(), modification_summary.as_str()),
            base_decoy_id: base_decoy.get_primary_key(),
            base_decoy: base_decoy.clone(),
            modification_ids: *Self::modifications_to_psql_array(modifications),
            modifications: modifications.clone(),
            weight: weight,
            modification_summary: modification_summary
        }
    }

    /// Creates a string of format "(modification_count|modification_accession|modification_description)..." for header. Will be prefixed with "ModSum=" in header.
    /// Sorted by "modification_accession|modification_description"
    fn create_modification_summary(modifications: &Vec<Modification>) -> String {
    // create a HashMap with key "mod_accession|mod_name" and value Vec<String> which contains the positions of the modification
        let mut modification_counts: HashMap<String, u8> = HashMap::new();
        for modification in modifications {
            let peff_notation_of_accession_and_name = format!("{}|{}", modification.get_accession(), modification.get_name());
            let counter = match modification_counts.entry(peff_notation_of_accession_and_name) {
                Vacant(entry) => entry.insert(0),
                Occupied(entry) => entry.into_mut()
            };
            *counter += 1;
        }
        let mut mod_res: String = String::new();
        let mut keys: Vec<String> = modification_counts.keys().map(|key| key.clone()).collect();
        keys.sort();    // sort keys, so summary modification is alway alphabetically sorted
        for peff_notation_of_accession_and_name in keys {
            if let Some(counter) = modification_counts.get(peff_notation_of_accession_and_name.as_str()) {
                mod_res.push_str(
                    format!(
                        "({}|{})",
                        counter,
                        peff_notation_of_accession_and_name
                    ).as_str()
                );
            }
        }
        return mod_res;
    }

    fn modifications_to_psql_array(modifications: &Vec<Modification>) -> Box<Array<i64>> {
        let modification_ids: Vec<i64> = modifications.iter().map(|modification| modification.get_primary_key()).collect();
        return Box::new(Array::from_vec(modification_ids, modifications.len() as i32));
    }

    fn modifications_from_psql_array(conn: &postgres::Connection, modification_ids: &Array<i64>) -> Result<Vec<Modification>, QueryError> {
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
        return Ok(modifications);
    }

    fn get_modification_ids(&self) -> &Array<i64> {
        return &self.modification_ids;
    }

    /// Creates Decoy from postgres::rows::Row. Created for `find_where_as_plain_decoys()`.
    ///
    /// # Arguments
    ///
    /// * `row` - row must contain [base_decoy.aa_sequence, base_decoy.header, modified_decoy.modification_summary, modified_decoy.weight] in this order
    fn plain_decoy_from_sql_row(row: &postgres::rows::Row) -> Decoy {
        return Decoy::new(
            format!("{} {}", row.get::<usize, String>(1), row.get::<usize, String>(2)).as_str(),
            &row.get::<usize, String>(0),
            &row.get(3)
        );
    }

    /// Works like find_where() but to save resources it does not create ModifiedDecoys with all Modification first but
    /// uses JOIN to get only attributes from BaseDecoys and ModifiedDecoys which are necessary for building a Decoy.
    /// Make sure that the number of $x used in `condition` are the same as elements in `values`.
    ///
    /// # Arguments
    ///
    /// * `conn` - Connection to Postgres
    /// * `condition` - WHERE-condition e.g. modified_decoy.weight BETWEEN $1 AND $2
    /// * `values` - e.g. [1000, 2000] for conditions example
    pub fn find_where_as_plain_decoys(conn: &postgres::Connection, conditions: &str, values: &[&postgres::types::ToSql]) -> Result<Vec<Decoy>, QueryError> {
        let select_query: String = format!(
            "SELECT {base_decoys_table}.aa_sequence, {base_decoys_table}.header, {modified_decoys_table}.modification_summary, {modified_decoys_table}.weight FROM {modified_decoys_table} INNER JOIN base_decoys ON {modified_decoys_table}.base_decoy_id={base_decoys_table}.id WHERE {conditions};",
            modified_decoys_table = Self::get_table_name(),
            base_decoys_table = BaseDecoy::get_table_name(),
            conditions = conditions
        );
        match conn.query(select_query.as_str(), values) {
            Ok(ref rows) => {
                let mut records: Vec<Decoy> = Vec::new();
                for row in rows {
                    records.push(Self::plain_decoy_from_sql_row(&row));
                }
                return Ok(records);
            },
            Err(err) => Err(handle_postgres_error(&err))
        }
    }
}

impl PeptideInterface for ModifiedDecoy {
    fn to_string(&self) -> String {
        return format!(
            "proteomic::modes::decoys::modified_decoy::ModifiedDecoy\n\theader => {}\n\taa_sequence => {}\n\tweight => {}",
            self.get_header(),
            self.get_aa_sequence(),
            mass::convert_mass_to_float(self.get_weight())
        );
    }

    fn get_header(&self) -> &str {
        return self.header.as_str();
    }

    fn get_aa_sequence(&self) -> String {
        return self.base_decoy.get_aa_sequence().to_owned();
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

        let modification_ids: Array<i64> = row.get::<usize, Array<i64>>(2);
        let modifications: Vec<Modification> = match Self::modifications_from_psql_array(&conn, &modification_ids) {
            Ok(modifications) => modifications,
            Err(err) => return Err(FromSqlRowError::InnerQueryError(err))
        };

        return Ok (
            Self {
                id: row.get(0),
                header: format!("{} ModSum={}", base_decoy.get_header(), row.get::<usize, String>(4)),
                base_decoy_id: base_decoy.get_primary_key(),
                base_decoy: base_decoy,
                modification_ids: modification_ids,
                modifications: modifications,
                weight: row.get(3),
                modification_summary: row.get(4)
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
        return "INSERT INTO modified_decoys (base_decoy_id, modification_ids, weight, modification_summary) VALUES ($1, $2, $3, $4) ON CONFLICT (base_decoy_id, modification_summary) DO NOTHING RETURNING id;";
    }

    fn create_attributes(&self) -> Box<Vec<&postgres::types::ToSql>>{
        return Box::new(vec![&self.base_decoy_id, &self.modification_ids, &self.weight, &self.modification_summary]);
    }

    fn update_query() -> &'static str{
        return "UPDATE modified_decoys SET base_decoy_id = $2, modification_ids = $3, weight = $4, modification_summary = $5 WHERE id = $1;";
    }

    fn update_attributes(&self) -> Box<Vec<&postgres::types::ToSql>>{
        return Box::new(vec![&self.id, &self.base_decoy_id, &self.modification_ids, &self.weight, &self.modification_summary]);
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
        return "SELECT id FROM modified_decoys WHERE base_decoy_id = $1 AND modification_summary = $2;"
    }

    fn exists_attributes(&self) -> Box<Vec<&postgres::types::ToSql>> {
        return Box::new(vec![&self.base_decoy_id, &self.modification_summary])
    }

    fn before_delete_hook(&self) -> Result<(), QueryError> {return Ok(());}
}

// PartialEq-implementation to use this type in a HashSet
impl PartialEq for ModifiedDecoy {
    fn eq(&self, other: &ModifiedDecoy) -> bool {
       return (self.get_aa_sequence() == other.get_aa_sequence())
        & (self.weight == other.get_weight())
        & (self.get_header() == other.get_header());
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
            modification_ids: self.modification_ids.clone(),
            modifications: self.modifications.clone(),
            weight: self.get_weight(),
            modification_summary: self.modification_summary.clone(),
            header: self.header.clone()
        }
    }
}