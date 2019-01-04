extern crate postgres;
extern crate postgres_array;

use std::hash::{Hash, Hasher};
use std::collections::HashMap;
use std::collections::hash_map::Entry::{Occupied, Vacant};

use self::postgres::Connection;
use self::postgres_array::Array;

use proteomic::models::decoys::decoy::Decoy;
use proteomic::models::decoys::base_decoy::BaseDecoy;
use proteomic::models::persistable::{handle_postgres_error, Persistable, QueryError, QueryOk, FromSqlRowError};
use proteomic::models::amino_acids::modification::Modification;
use proteomic::models::mass;
use proteomic::utility::database_connection::DatabaseConnection;

pub struct ModifiedDecoy {
    id: i64,                                        // SERIAL
    base_decoy: BaseDecoy,                          // BIGINT REFERENCE
    c_terminus_modification: Option<Modification>,  // BIGINT REFERENCE
    n_terminus_modification: Option<Modification>,  // BIGINT REFERENCE
    modifications: Vec<Option<Modification>>,       // INTEGER ARRAY
    weight: i64                                     // BIGINT
}

impl ModifiedDecoy {
    pub fn new(base_decoy: &BaseDecoy, c_terminus_modification: Option<Modification>, n_terminus_modification: Option<Modification>, modifications: &Vec<Option<Modification>>, weight: i64) -> ModifiedDecoy {
        return Self {
            id: 0,
            base_decoy: base_decoy.clone(),
            c_terminus_modification: c_terminus_modification,
            n_terminus_modification: n_terminus_modification,
            modifications: modifications.clone(),
            weight: weight
        }
    }

    fn modifications_to_psql_array(&self) -> Box<Array<Option<i64>>> {
        let mut modification_ids: Vec<Option<i64>> = Vec::new();
        for modification_option in self.modifications.iter() {
            match modification_option {
                Some(modification) => modification_ids.push(Some(modification.get_primary_key())),
                None => modification_ids.push(None)
            }
        }
        return Box::new(Array::from_vec(modification_ids, self.base_decoy.get_length()));
    }

    pub fn get_c_terminus_modification(&self) -> &Option<Modification> {
        return &self.c_terminus_modification;
    }

    pub fn get_n_terminus_modification(&self) -> &Option<Modification> {
        return &self.n_terminus_modification
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
        // create a HashMap with key "mod_accession|mod_name" and value Vec<String> which contains the positions of the modification
        let mut accession_header_position_map: HashMap<String, Vec<String>> = HashMap::new();
        for idx in 0..(self.modifications.len()) {
            match &self.modifications[idx] {
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
        let mut mod_res_unimod: String = String::new();
        for (accession_and_header, positions) in &accession_header_position_map {
            mod_res_unimod.push_str(
                format!(
                    "({}|{})",
                    positions.join(","),
                    accession_and_header
                ).as_str()
            );
        }
        return format!(
            "{} ModResUniMod={}",
            self.base_decoy.get_header(),
            mod_res_unimod
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


impl Persistable<ModifiedDecoy, i64, (i64, Option<i64>, Option<i64>, &Array<Option<i64>>, i64)> for ModifiedDecoy {
    fn from_sql_row(row: &postgres::rows::Row) -> Result<Self, FromSqlRowError> {
        let conn = DatabaseConnection::get_database_connection();
        let base_decoy: BaseDecoy = match BaseDecoy::find(&conn, &row.get::<usize, i64>(1)) {
            Ok(base_decoy) => base_decoy,
            Err(query_error) => match query_error {
                QueryError::NoMatch => return Err(FromSqlRowError::AssociatedRecordNotFound(format!("BaseDecoy(id: {})", row.get::<usize, i64>(1)))),
                _ => return Err(FromSqlRowError::InnerQueryError(query_error))                          // else a sql error occures
            }
        };
        return Ok (
            Self {
                modifications: vec![None; base_decoy.get_length() as usize],
                id: row.get(0),
                base_decoy: base_decoy,
                c_terminus_modification: match Modification::find(&conn, &row.get::<usize, i64>(2)) {
                    Ok(modification) => Some(modification),
                    Err(error) => match error {
                        QueryError::NoMatch => None,                                // if error is NOHIT there is simply no record
                        _ => return Err(FromSqlRowError::InnerQueryError(error))  // else a sql error occures
                    }
                },
                n_terminus_modification: match Modification::find(&conn, &row.get::<usize, i64>(2)) {
                    Ok(modification) => Some(modification),
                    Err(error) => match error {
                        QueryError::NoMatch => None,                                    // if error is NOHIT there is simply no record
                        _ => return Err(FromSqlRowError::InnerQueryError(error))      // else a sql error occures
                    }
                },
                weight: row.get(4)
            }
        );
    }

    fn get_primary_key(&self) -> i64 {
        return self.id;
    }

    fn find(conn: &Connection, primary_key: &i64) -> Result<Self, QueryError> {
        match conn.query("SELECT * FROM modified_decoys WHERE id = $1 LIMIT 1", &[primary_key]) {
            Ok(ref rows) if rows.len() > 0 => match Self::from_sql_row(&rows.get(0)) {
                Ok(record) => Ok(record),
                Err(err) => match err {
                    FromSqlRowError::InnerQueryError(from_sql_err) => Err(from_sql_err),
                    FromSqlRowError::AssociatedRecordNotFound(from_sql_err) => Err(QueryError::AssociatedRecordNotFound(from_sql_err.to_string()))
                }
            },
            Ok(_rows) => Err(QueryError::NoMatch),
            Err(err) => Err(handle_postgres_error(&err))
        }
    }

    fn find_by_unique_identifier(conn: &Connection, unique_identifier: &(i64, Option<i64>, Option<i64>, &Array<Option<i64>>, i64)) -> Result<Self, QueryError> {
        match conn.query(
            "SELECT * FROM modified_decoys WHERE base_decoy_id = $1 AND c_terminus_modification_id = $2 AND c_terminus_modification_id = $3 AND modification_ids = $4 AND weight = $5 LIMIT 1",
            &[&unique_identifier.0, &unique_identifier.1, &unique_identifier.2, &unique_identifier.3, &unique_identifier.4]
        ) {
            Ok(ref rows) if rows.len() > 0 => match Self::from_sql_row(&rows.get(0)) {
                Ok(record) => Ok(record),
                Err(err) => match err {
                    FromSqlRowError::InnerQueryError(from_sql_err) => Err(from_sql_err),
                    FromSqlRowError::AssociatedRecordNotFound(from_sql_err) => Err(QueryError::AssociatedRecordNotFound(from_sql_err.to_string()))
                }
            },
            Ok(_rows) => Err(QueryError::NoMatch),
            Err(err) => Err(handle_postgres_error(&err))
        }
    }


    fn create(&mut self, conn: &postgres::Connection) -> Result<QueryOk, QueryError> {
        // Option<i64> is nullable for sql see: https://github.com/sfackler/rust-postgres#type-correspondence
        let c_terminus_modification_id: Option<i64> = match &self.c_terminus_modification {
            Some(modification) => Some(modification.get_primary_key()),
            None => None
        };
        // Option<i64> is nullable for sql see: https://github.com/sfackler/rust-postgres#type-correspondence
        let n_terminus_modification_id: Option<i64> = match &self.c_terminus_modification {
            Some(modification) => Some(modification.get_primary_key()),
            None => None
        };
        let modification_ids_boxed = self.modifications_to_psql_array();
        match conn.query(
            Self::get_insert_query(),
            &[&self.base_decoy.get_primary_key(), &c_terminus_modification_id, &n_terminus_modification_id, modification_ids_boxed.as_ref(), &self.weight]
        ) {
            Ok(ref rows) if rows.len() > 0 => {
                self.id = rows.get(0).get(0);
                return Ok(QueryOk::Created);
            },
            Ok(_rows) => {
                // zero rows means there are a conflict on update, so the decoys exists already
                match Self::find_by_unique_identifier(conn, &(self.base_decoy.get_primary_key(), c_terminus_modification_id, n_terminus_modification_id, modification_ids_boxed.as_ref(), self.weight)) {
                    Ok(modified_decoy) => {
                        self.id = modified_decoy.get_primary_key();
                        return Ok(QueryOk::AlreadyExists);
                    },
                    Err(err) => Err(err)
                }
            }
            Err(err) => Err(handle_postgres_error(&err))
        }
    }

    fn update(&mut self, conn: &postgres::Connection) -> Result<QueryOk, QueryError> {
        if !self.is_persisted() {
            return Err(QueryError::RecordIsNotPersisted);
        }
        // Option<i64> is nullable for sql see: https://github.com/sfackler/rust-postgres#type-correspondence
        let c_terminus_modification_id: Option<i64> = match &self.c_terminus_modification {
            Some(modification) => Some(modification.get_primary_key()),
            None => None
        };
        // Option<i64> is nullable for sql see: https://github.com/sfackler/rust-postgres#type-correspondence
        let n_terminus_modification_id: Option<i64> = match &self.c_terminus_modification {
            Some(modification) => Some(modification.get_primary_key()),
            None => None
        };
        let modification_ids_boxed = self.modifications_to_psql_array();
        match conn.query(
            Self::get_update_query(),
            &[&self.id, &self.base_decoy.get_primary_key(), &c_terminus_modification_id, &n_terminus_modification_id, modification_ids_boxed.as_ref(), &self.weight]
        ) {
            Ok(ref rows) if rows.len() > 0 => Ok(QueryOk::Updated),
            Ok(_rows) => Err(QueryError::NoReturn),
            Err(err) => Err(handle_postgres_error(&err))
        }
    }

    fn save(&mut self, conn: &postgres::Connection) -> Result<QueryOk, QueryError> {
        if self.is_persisted() {
            return self.update(conn);
        } else {
            return self.create(conn);
        }
    }

    fn delete(&mut self, conn: &postgres::Connection) -> Result<QueryOk, QueryError> {
        if !self.is_persisted() {
            return Err(QueryError::RecordIsNotPersisted);
        }
        match conn.execute("DELETE FROM modified_decoys WHERE id = $1;", &[&self.id]) {
            Ok(_) => {
                self.id = 0;
                return Ok(QueryOk::Deleted);
            },
            Err(err) => Err(handle_postgres_error(&err))
        }
    }

    fn delete_all(conn: &postgres::Connection) -> Result<QueryOk, QueryError> {
        match conn.execute("DELETE FROM modified_decoys WHERE id IS NOT NULL;", &[]) {
            Ok(_) => Ok(QueryOk::Deleted),
            Err(err) => Err(handle_postgres_error(&err))
        }
    }


    fn get_table_name() -> &'static str {
        return "modified_decoys";
    }

    fn get_select_primary_key_by_unique_identifier_query() -> &'static str {
        return "SELECT * FROM modified_decoys WHERE base_decoy_id = $1 AND c_terminus_modification_id = $2 AND c_terminus_modification_id = $3 AND modification_ids = $4 AND weight = $5 LIMIT 1";
    }

    fn get_insert_query() -> &'static str {
        return "INSERT INTO modified_decoys (base_decoy_id, c_terminus_modification_id, n_terminus_modification_id, modification_ids, weight) VALUES ($1, $2, $3, $4, $5) ON CONFLICT (base_decoy_id, c_terminus_modification_id, n_terminus_modification_id, modification_ids, weight) DO NOTHING RETURNING id";
    }

    fn get_update_query() -> &'static str {
        return "UPDATE modified_decoys SET base_decoy_id = $2, c_terminus_modification_id = $3, n_terminus_modification_id = $4, modification_ids = $5, weight = $6 WHERE id = $1";
    }

    fn exec_select_primary_key_by_unique_identifier_statement(&mut self, prepared_statement: &postgres::stmt::Statement) -> Result<QueryOk, QueryError> {
        // Option<i64> is nullable for sql see: https://github.com/sfackler/rust-postgres#type-correspondence
        let c_terminus_modification_id: Option<i64> = match &self.c_terminus_modification {
            Some(modification) => Some(modification.get_primary_key()),
            None => None
        };
        // Option<i64> is nullable for sql see: https://github.com/sfackler/rust-postgres#type-correspondence
        let n_terminus_modification_id: Option<i64> = match &self.c_terminus_modification {
            Some(modification) => Some(modification.get_primary_key()),
            None => None
        };
        match prepared_statement.query(&[&self.base_decoy.get_primary_key(), &c_terminus_modification_id, &n_terminus_modification_id, self.modifications_to_psql_array().as_ref(), &self.weight]) {
            Ok(ref rows) if rows.len() > 0 => {
                self.id = rows.get(0).get(0);
                return Ok(QueryOk::Selected);
            },
            Ok(_rows) => Err(QueryError::NoMatch),
            Err(err) => Err(handle_postgres_error(&err))
        }
    }

    fn exec_insert_statement(&mut self, prepared_statement: &postgres::stmt::Statement) -> Result<QueryOk, QueryError> {
        // Option<i64> is nullable for sql see: https://github.com/sfackler/rust-postgres#type-correspondence
        let c_terminus_modification_id: Option<i64> = match &self.c_terminus_modification {
            Some(modification) => Some(modification.get_primary_key()),
            None => None
        };
        // Option<i64> is nullable for sql see: https://github.com/sfackler/rust-postgres#type-correspondence
        let n_terminus_modification_id: Option<i64> = match &self.c_terminus_modification {
            Some(modification) => Some(modification.get_primary_key()),
            None => None
        };
        match prepared_statement.query(&[&self.base_decoy.get_primary_key(), &c_terminus_modification_id, &n_terminus_modification_id, self.modifications_to_psql_array().as_ref(), &self.weight]) {
            Ok(ref rows) if rows.len() > 0 => {
                self.id = rows.get(0).get(0);
                return Ok(QueryOk::Created);
            },
            Ok(_rows) => Err(QueryError::NoReturn),
            Err(err) => Err(handle_postgres_error(&err))
        }
    }

    fn exec_update_statement(&mut self, prepared_statement: &postgres::stmt::Statement) -> Result<QueryOk, QueryError> {
        // Option<i64> is nullable for sql see: https://github.com/sfackler/rust-postgres#type-correspondence
        let c_terminus_modification_id: Option<i64> = match &self.c_terminus_modification {
            Some(modification) => Some(modification.get_primary_key()),
            None => None
        };
        // Option<i64> is nullable for sql see: https://github.com/sfackler/rust-postgres#type-correspondence
        let n_terminus_modification_id: Option<i64> = match &self.c_terminus_modification {
            Some(modification) => Some(modification.get_primary_key()),
            None => None
        };
        match prepared_statement.query(&[&self.id, &self.base_decoy.get_primary_key(), &c_terminus_modification_id, &n_terminus_modification_id, self.modifications_to_psql_array().as_ref(), &self.weight]) {
            Ok(ref rows) if rows.len() > 0 => Ok(QueryOk::Updated),
            Ok(_rows) => Err(QueryError::NoReturn),
            Err(err) => Err(handle_postgres_error(&err))
        }
    }


    fn is_persisted(&self) -> bool {
        return self.id > 0;
    }
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
            base_decoy: self.base_decoy.clone(),
            c_terminus_modification: self.c_terminus_modification.clone(),
            n_terminus_modification: self.n_terminus_modification.clone(),
            modifications: self.modifications.clone(),
            weight: self.get_weight()
        }
    }
}