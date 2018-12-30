extern crate postgres;
extern crate postgres_array;

use std::hash::{Hash, Hasher};

use self::postgres::Connection;
use self::postgres_array::Array;

use proteomic::models::decoys::decoy::Decoy;
use proteomic::models::decoys::base_decoy::BaseDecoy;
use proteomic::models::persistable::Persistable;
use proteomic::models::amino_acids::modification::{Modification, ModificationPosition};
use proteomic::models::mass;

pub struct ModifiedDecoy {
    id: i64,                                        // SERIAL
    base_decoy: BaseDecoy,                          // BIGINT REFERENCE
    c_terminus_modification: Option<Modification>,  // BIGINT REFERENCE
    n_terminus_modification: Option<Modification>,  // BIGINT REFERENCE
    modifications: Vec<Option<Modification>>,       // INTEGER ARRAY
    weight: i64                                     // BIGINT
}

impl ModifiedDecoy {
    pub fn new(base_decoy: BaseDecoy) -> ModifiedDecoy {
        return Self {
            modifications: vec![None; base_decoy.get_length() as usize],
            weight: base_decoy.get_weight(),
            id: -1,
            base_decoy: base_decoy,
            c_terminus_modification: None,
            n_terminus_modification: None
        }
    }

    pub fn new_from_psql_rows(conn: &Connection, rows: &postgres::rows::Rows) -> Result<Self, String> {
        let base_decoy: BaseDecoy = match BaseDecoy::find(conn, &rows.get(0).get::<usize, i64>(1)) {
            Ok(base_decoy) => base_decoy,
            Err(error) => match error.as_str() {
                "NOHIT" => return Err(format!("could not found the associated BaseDecoy with id = {}", rows.get(0).get::<usize, i64>(1))),
                _ => return Err(error)                          // else a sql error occures
            }
        };
        return Ok (
            Self {
                modifications: vec![None; base_decoy.get_length() as usize],
                id: rows.get(0).get(0),
                base_decoy: base_decoy,
                c_terminus_modification: match Modification::find(conn, &rows.get(0).get::<usize, i64>(2)) {
                    Ok(modification) => Some(modification),
                    Err(error) => match error.as_str() {
                        "NOHIT" => None,                            // if error is NOHIT there is simply no record
                        _ => return Err(error)                      // else a sql error occures
                    }
                },
                n_terminus_modification: match Modification::find(conn, &rows.get(0).get::<usize, i64>(2)) {
                    Ok(modification) => Some(modification),
                    Err(error) => match error.as_str() {
                        "NOHIT" => None,                            // if error is NOHIT there is simply no record
                        _ => return Err(error)                      // else a sql error occures
                    }
                },
                weight: rows.get(0).get(4)
            }
        );
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

    fn set_c_terminus_modification_to_none(&mut self) {
        match &self.c_terminus_modification {
            Some(modification) => self.weight -= modification.get_mono_mass(),
            None => ()
        };
        self.c_terminus_modification = None;
    }

    fn set_n_terminus_modification_to_none(&mut self) {
        match &self.n_terminus_modification {
            Some(modification) => self.weight -= modification.get_mono_mass(),
            None => ()
        };
        self.n_terminus_modification = None;
    }

    fn set_modification_to_none_at(&mut self, idx: usize) {
        if idx < self.modifications.len() {
            match &self.modifications.get(idx) {
                Some(modification_option) => {
                    match modification_option {
                        Some(modification) => self.weight -= modification.get_mono_mass(),
                        None => ()
                    }
                }
                None => ()
            };
            self.modifications[idx] = None;
        }
    }

    pub fn set_c_terminus_modification(&mut self, c_terminus_modification: Option<Modification>) -> bool {
        match c_terminus_modification {
            Some(modification) => {
                if (modification.get_position() == ModificationPosition::CTerminus) | (modification.get_position() == ModificationPosition::Anywhere) {
                    self.set_c_terminus_modification_to_none();
                    self.weight += modification.get_mono_mass();
                    self.c_terminus_modification = Some(modification);
                    return true;
                }
            }
            None => ()
        }
        return false;
    }

    pub fn set_n_terminus_modification(&mut self, n_terminus_modification: Option<Modification>) -> bool {
        match n_terminus_modification {
            Some(modification) => {
                if (modification.get_position() == ModificationPosition::NTerminus) | (modification.get_position() == ModificationPosition::Anywhere) {
                    self.set_n_terminus_modification_to_none();
                    self.weight += modification.get_mono_mass();
                    self.n_terminus_modification = Some(modification);
                    return true;
                }
            }
            None => ()
        }
        return false;
    }

    pub fn set_modification_at(&mut self, idx: usize, modification_option: Option<Modification>) -> bool {
        match modification_option {
            Some(modification) => {
                if (idx < self.modifications.len()) & (modification.get_position() == ModificationPosition::Anywhere) {
                    self.set_modification_to_none_at(idx);
                    self.weight += modification.get_mono_mass();
                    self.modifications[idx] = Some(modification);
                    return true;
                }
            }
            None => ()
        }
        return false;
    }
}

impl Decoy for ModifiedDecoy {
    fn to_string(&self) -> String {
        return format!(
            "proteomic::modes::decoys::base_decoy::ModifiedDecoy\n\taa_sequence => {}\n\tweight => {}",
            self.get_aa_sequence(),
            mass::convert_mass_to_float(self.get_weight())
        );
    }

    fn get_header(&self) -> String {
        // add modifications to header here and return
        return self.base_decoy.get_header();
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
}


impl Persistable<ModifiedDecoy, i64, (i64, Option<i64>, Option<i64>, &Array<Option<i64>>, i64)> for ModifiedDecoy {
    fn get_primary_key(&self) -> i64 {
        return self.id;
    }

    fn find(conn: &Connection, primary_key: &i64) -> Result<Self, String> {
        match conn.query("SELECT * FROM base_decoys WHERE id = $1 LIMIT 1", &[primary_key]) {
            Ok(rows) =>{
                if rows.len() > 0 {
                    match Self::new_from_psql_rows(conn, &rows) {
                        Ok(modified_decoy) => Ok(modified_decoy),
                        Err(err) => Err(err)
                    }
                } else {
                    Err("NOHIT".to_owned())
                }
            },
            Err(err) => Err(err.code().unwrap().code().to_owned())
        }
    }

    fn find_by_unique_identifier(conn: &Connection, unique_identifier: &(i64, Option<i64>, Option<i64>, &Array<Option<i64>>, i64)) -> Result<Self, String> {
        match conn.query(
            "SELECT * FROM base_decoys WHERE base_decoy_id = $1 AND c_terminus_modification_id = $2 AND c_terminus_modification_id = $3 AND modification_ids = $4 AND weight = $5 LIMIT 1",
            &[&unique_identifier.0, &unique_identifier.1, &unique_identifier.2, &unique_identifier.3, &unique_identifier.4]
        ) {
            Ok(ref rows) if rows.len() > 0 => {
                match Self::new_from_psql_rows(conn, &rows) {
                    Ok(modified_decoy) => Ok(modified_decoy),
                    Err(err) => Err(err)
                }
            },
            Ok(_rows) => Err("NOHIT".to_owned()),
            Err(err) => Err(err.code().unwrap().code().to_owned())
        }
    }


    fn create(&mut self, conn: &postgres::Connection) -> Result<(), String> {
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
                self.id =  rows.get(0).get(0);
                return Ok(());
            },
            Ok(_rows) => {
                // zero rows means there are a conflict on update, so the decoys exists already
                match Self::find_by_unique_identifier(conn, &(self.base_decoy.get_primary_key(), c_terminus_modification_id, n_terminus_modification_id, modification_ids_boxed.as_ref(), self.weight)) {
                    Ok(modified_decoy) => {
                        self.id = modified_decoy.get_primary_key();
                        return Ok(());
                    },
                    Err(err) => Err(format!("cannot insert nor find ModifiedDecoy\n\toriginal error: {}", err))
                }
            }
            Err(err) => Err(err.code().unwrap().code().to_owned())
        }
    }

    fn update(&mut self, conn: &postgres::Connection) -> Result<(), String> {
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
            Ok(ref rows) if rows.len() > 0 => Ok(()),
            Ok(_rows) => Err("NORET".to_owned()),
            Err(err) => Err(err.code().unwrap().code().to_owned())
        }
    }

    fn save(&mut self, conn: &postgres::Connection) -> Result<(), String> {
        if self.is_persisted() {
            return self.update(conn);
        } else {
            return self.create(conn);
        }
    }


    fn get_select_primary_key_by_unique_identifier_query() -> &'static str {
        return "SELECT * FROM base_decoys WHERE base_decoy_id = $1 AND c_terminus_modification_id = $2 AND c_terminus_modification_id = $3 AND modification_ids = $4 AND weight = $5 LIMIT 1";
    }

    fn get_insert_query() -> &'static str {
        return "INSERT INTO base_decoys (base_decoy_id, c_terminus_modification_id, n_terminus_modification_id, modification_ids, weight) VALUES ($1, $2, $3, $4) ON CONFLICT (base_decoy_id, c_terminus_modification_id, n_terminus_modification_id, modification_ids, weight) DO NOTHING RETURNING id";
    }

    fn get_update_query() -> &'static str {
        return "UPDATE base_decoys SET base_decoy_id = $2, c_terminus_modification_id = $3, n_terminus_modification_id = $4, modification_ids = $5, weight = $6 WHERE id = $1";
    }

    fn exec_select_primary_key_by_unique_identifier_statement(&mut self, prepared_statement: &postgres::stmt::Statement) -> Result<(), String> {
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
                return Ok(());
            },
            Ok(_rows) => Err("NOHIT".to_owned()),
            Err(err) => Err(err.code().unwrap().code().to_owned())
        }
    }

    fn exec_insert_statement(&mut self, prepared_statement: &postgres::stmt::Statement) -> Result<(), String> {
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
                return Ok(());
            },
            Ok(_rows) => Err("NORET".to_owned()),
            Err(err) => Err(err.code().unwrap().code().to_owned())
        }
    }

    fn exec_update_statement(&mut self, prepared_statement: &postgres::stmt::Statement) -> Result<(), String> {
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
            Ok(ref rows) if rows.len() > 0 => Ok(()),
            Ok(_rows) => Err("NORET".to_owned()),
            Err(err) => Err(err.code().unwrap().code().to_owned())
        }
    }


    fn is_persisted(&self) -> bool {
        return self.id > 0;
    }

    fn get_count(conn: &postgres::Connection) -> i64 {
        return match conn.query("SELECT cast(count(id) AS BIGINT) FROM decoys", &[]) {
            Ok(ref rows) if rows.len() > 0 => rows.get(0).get::<usize, i64>(0),
            _ => -1
        };
    }

}


// PartialEq-implementation to use this type in a HashSet
impl PartialEq for ModifiedDecoy {
    fn eq(&self, other: &ModifiedDecoy) -> bool {
       return (self.get_aa_sequence() == *other.get_aa_sequence()) & (self.weight == other.get_weight());
    }
}

// Eq-implementation to use this type in a HashSet
impl Eq for ModifiedDecoy {}

// Hash-implementation to use this type in a HashSet
impl Hash for ModifiedDecoy {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.get_aa_sequence().hash(state);
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