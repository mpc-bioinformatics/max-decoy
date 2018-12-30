extern crate postgres;

use std::hash::{Hash, Hasher};

use self::postgres::Connection;

use proteomic::models::decoys::decoy::Decoy;
use proteomic::models::persistable::Persistable;
use proteomic::models::peptide::Peptide;
use proteomic::models::mass;

pub struct BaseDecoy {
    id: i64,                            // BIGSERIAL
    header: String,                     // TEXT
    aa_sequence: String,                // VARCHAR(60)
    length: i32,                        // INTEGER
    weight: i64                         // BIGINT
}

impl BaseDecoy {
    pub fn new(header: &str, aa_sequence: &str, weight: i64) -> BaseDecoy {
        return Self {
            id: -1,
            header: header.to_owned(),
            length: aa_sequence.len() as i32,
            aa_sequence: aa_sequence.to_owned(),
            weight: weight
        }
    }
}

impl Decoy for BaseDecoy {
    fn to_string(&self) -> String {
        return format!(
            "proteomic::modes::decoys::base_decoy::BaseDecoy\n\taa_sequence => {}\n\tweight => {}", 
            self.aa_sequence, 
            mass::convert_mass_to_float(self.weight)
        );
    }

    fn get_header(&self) -> String {
        return self.header.clone();
    }

    fn get_aa_sequence(&self) -> String {
        return self.aa_sequence.clone();
    }

    fn get_weight(&self) -> i64 {
        return self.weight;
    }

    fn get_length(&self) -> i32 {
        return self.aa_sequence.len() as i32;
    }

    fn get_c_terminus_amino_acid(&self) -> char {
        match self.aa_sequence.chars().last() {
            Some(amino_acids_one_letter_code) => amino_acids_one_letter_code,
            None => '_'
        }
    }

    fn get_n_terminus_amino_acid(&self) -> char {
        match self.aa_sequence.chars().next() {
            Some(amino_acids_one_letter_code) => amino_acids_one_letter_code,
            None => '_'
        }
    }

    fn get_amino_acid_at(&self, idx: usize) -> char {
        match self.aa_sequence.chars().nth(idx) {
            Some(amino_acids_one_letter_code) => amino_acids_one_letter_code,
            None => '_'
        }
    }
}


impl Persistable<BaseDecoy, i64, String> for BaseDecoy {
    fn from_sql_row(row: &postgres::rows::Row) -> Result<Self, String> {
        return Ok(
            Self {
                id: row.get(0),
                header: row.get(1),
                aa_sequence: row.get(2),
                length: row.get(3),
                weight: row.get(4)
            }
        )
    }


    fn get_primary_key(&self) -> i64 {
        return self.id;
    }

    fn find(conn: &Connection, primary_key: &i64) -> Result<Self, String> {
        match conn.query("SELECT * FROM base_decoys WHERE id = $1 LIMIT 1", &[primary_key]) {
            Ok(ref rows) if rows.len() > 0 => Self::from_sql_row(&rows.get(0)),
            Ok(_rows) => Err("NOHIT".to_owned()),
            Err(err) => Err(err.code().unwrap().code().to_owned())
        }
    }

    fn find_by_unique_identifier(conn: &Connection, unique_identifier: &String) -> Result<Self, String> {
        let generalized_aa_sequence: String = Peptide::gerneralize_aa_sequence(unique_identifier);
        match conn.query(
            "SELECT * FROM base_decoys WHERE aa_sequence = $1 LIMIT 1",
            &[&generalized_aa_sequence]
        ) {
            Ok(ref rows) if rows.len() > 0 => Self::from_sql_row(&rows.get(0)),
            Ok(_rows) => Err("NOHIT".to_owned()),
            Err(err) => Err(err.code().unwrap().code().to_owned())
        }
    }


    fn create(&mut self, conn: &postgres::Connection) -> Result<(), String> {
        match conn.query(
            Self::get_insert_query(),
            &[&self.header, &self.aa_sequence, &self.weight, &self.length]
        ) {
            Ok(ref rows) if rows.len() > 0 => {
                self.id =  rows.get(0).get(0);
                return Ok(());
            },
            Ok(_rows) => {
                // zero rows means there are a conflict on update, so the decoys exists already
                match Self::find_by_unique_identifier(conn, &self.aa_sequence) {
                    Ok(decoy) => {
                        self.id = decoy.get_primary_key();
                        return Ok(());
                    },
                    Err(err) => Err(format!("cannot insert nor find decoy '{}'\n\toriginal error: {}", self.aa_sequence, err))
                }
            }
            Err(err) => Err(err.code().unwrap().code().to_owned())
        }
    }

    fn update(&mut self, conn: &postgres::Connection) -> Result<(), String> {
        match conn.query(
            Self::get_update_query(),
            &[&self.id, &self.header, &self.aa_sequence, &self.weight, &self.length]
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

    fn delete(&mut self, conn: &postgres::Connection) -> Result<(), String> {
        if !self.is_persisted() {
            return Err("BaseDecoy is not persisted".to_owned());
        }
        match conn.execute("DELETE FROM base_decoys WHERE id = $1;", &[&self.id]) {
            Ok(_) => {
                self.id = 0;
                return Ok(());
            },
            Err(err) => Err(format!("could not delete BaseDecoy from database; postgresql error is: {}", err))
        }
    }

    fn delete_all(conn: &postgres::Connection) -> Result<(), String> {
        match conn.execute("DELETE FROM base_decoys WHERE id IS NOT NULL;", &[]) {
            Ok(_) => Ok(()),
            Err(err) => Err(format!("could not delete BaseDecoys from database; postgresql error is: {}", err))
        }
    }

    fn select_where(conn: &postgres::Connection, conditions: &str, values: &[&postgres::types::ToSql]) -> Result<Vec<Self>, String> {
        let where_statement: String = format!("SELECT * FROM base_decoys WHERE {};", conditions);
        match conn.query(where_statement.as_str(), values) {
            Ok(ref rows) => {
                let records: Vec<Self> = Vec::new();
                for row in rows {
                    match Self::from_sql_row(&row) {
                        Ok(record) => records.push(record),
                        Err(err) => return Err(err)
                    }
                    
                }
                return Ok(records);
            },
            Err(err) => Err(format!("could not gether record from database; postgresql error is: {}", err))
        }
    }


    fn get_table_name() -> &'static str {
        return "base_decoys";
    }

    fn get_select_primary_key_by_unique_identifier_query() -> &'static str {
        return "SELECT id FROM base_decoys WHERE aa_sequence = $1 LIMIT 1";
    }

    fn get_insert_query() -> &'static str {
        return "INSERT INTO base_decoys (header, aa_sequence, weight, length) VALUES ($1, $2, $3, $4) ON CONFLICT (aa_sequence, weight) DO NOTHING RETURNING id";
    }

    fn get_update_query() -> &'static str {
        return "UPDATE base_decoys SET header = $2, aa_sequence = $3, weight = $4, length = $5 WHERE id = $1";
    }

    fn exec_select_primary_key_by_unique_identifier_statement(&mut self, prepared_statement: &postgres::stmt::Statement) -> Result<(), String> {
        match prepared_statement.query(&[&self.aa_sequence]) {
            Ok(ref rows) if rows.len() > 0 => {
                self.id = rows.get(0).get(0);
                return Ok(());
            },
            Ok(_rows) => Err("NOHIT".to_owned()),
            Err(err) => Err(err.code().unwrap().code().to_owned())
        }
    }

    fn exec_insert_statement(&mut self, prepared_statement: &postgres::stmt::Statement) -> Result<(), String> {
        match prepared_statement.query(&[&self.aa_sequence, &self.weight, &self.length]) {
            Ok(ref rows) if rows.len() > 0 => {
                self.id = rows.get(0).get(0);
                return Ok(());
            },
            Ok(_rows) => Err("NORET".to_owned()),
            Err(err) => Err(err.code().unwrap().code().to_owned())
        }
    }

    fn exec_update_statement(&mut self, prepared_statement: &postgres::stmt::Statement) -> Result<(), String> {
        match prepared_statement.query(&[&self.id, &self.aa_sequence, &self.weight, &self.length]) {
            Ok(ref rows) if rows.len() > 0 => Ok(()),
            Ok(_rows) => Err("NORET".to_owned()),
            Err(err) => Err(err.code().unwrap().code().to_owned())
        }
    }


    fn is_persisted(&self) -> bool {
        return self.id > 0;
    }
}


// PartialEq-implementation to use this type in a HashSet
impl PartialEq for BaseDecoy {
    fn eq(&self, other: &BaseDecoy) -> bool {
       return (self.aa_sequence == *other.get_aa_sequence()) & (self.weight == other.get_weight());
    }
}

// Eq-implementation to use this type in a HashSet
impl Eq for BaseDecoy {}

// Hash-implementation to use this type in a HashSet
impl Hash for BaseDecoy {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.aa_sequence.hash(state);
    }
}

impl Clone for BaseDecoy {
    fn clone(&self) -> Self {
        return Self{
            id: self.id,
            header: self.header.clone(),
            aa_sequence: String::from(self.aa_sequence.as_str()),
            length: self.length,
            weight: self.weight
        }
    }
}