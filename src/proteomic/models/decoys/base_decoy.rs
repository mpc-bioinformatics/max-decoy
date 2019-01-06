extern crate postgres;

use std::hash::{Hash, Hasher};

use proteomic::models::decoys::decoy::{Decoy, PlainDecoy};
use proteomic::models::persistable::{handle_postgres_error, Persistable, QueryError, FromSqlRowError};
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

    // row must contain [header, aa_sequence, weight] in this order
    fn prepared_decoy_from_sql_row(row: &postgres::rows::Row) -> PlainDecoy {
        return PlainDecoy::new(
            format!("{} {}", row.get::<usize, String>(1), row.get::<usize, String>(2)).as_str(),
            &row.get::<usize, String>(0), 
            &row.get(3)
        );
    }

    // this function works like find_where() but to save time it will return PlainDecoy instead of loading all modifications first
    pub fn find_where_as_prepared_decoys(conn: &postgres::Connection, conditions: &str, values: &[&postgres::types::ToSql]) -> Result<Vec<PlainDecoy>, QueryError> {
        let select_query: String = format!("SELECT header, aa_sequence, weight FROM {} WHERE {};", Self::get_table_name(), conditions);
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

impl Decoy for BaseDecoy {
    fn to_string(&self) -> String {
        return format!(
            "proteomic::modes::decoys::base_decoy::BaseDecoy\n\theader => {}\n\taa_sequence => {}\n\tweight => {}",
            self.header,
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
    fn from_sql_row(row: &postgres::rows::Row) -> Result<Self, FromSqlRowError> {
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
        return "base_decoys";
    }

    fn is_persisted(&self) -> bool {
        return self.id > 0;
    }

    fn find_query() -> &'static str {
        return "SELECT * FROM base_decoys WHERE id = $1 LIMIT 1;";
    }

    fn create_query() -> &'static str {
        return "INSERT INTO base_decoys (header, aa_sequence, weight, length) VALUES ($1, $2, $3, $4) ON CONFLICT (aa_sequence, weight) DO NOTHING RETURNING id;";
    }

    fn create_attributes(&self) -> Box<Vec<&postgres::types::ToSql>>{
        return Box::new(vec![&self.header, &self.aa_sequence, &self.weight, &self.length]);
    }

    fn update_query() -> &'static str{
        return "UPDATE base_decoys SET header = $2, aa_sequence = $3, weight = $4, length = $5 WHERE id = $1;";
    }

    fn update_attributes(&self) -> Box<Vec<&postgres::types::ToSql>>{
        return Box::new(vec![&self.id, &self.header, &self.aa_sequence, &self.weight, &self.length])
    }

    fn delete_query() -> &'static str {
        return "DELETE FROM base_decoys WHERE id = $1;";
    }

    fn delete_attributes(&self) -> Box<Vec<&postgres::types::ToSql>> {
        return Box::new(vec![&self.id]);
    }

    fn delete_all_query() -> &'static str {
        return "DELETE FROM base_decoys WHERE id IS NOT NULL;";
    }

    fn exists_query() -> &'static str {
        return "SELECT id FROM base_decoys WHERE aa_sequence = $1 LIMIT 1;";
    }

    fn exists_attributes(&self) -> Box<Vec<&postgres::types::ToSql>> {
        return Box::new(vec![&self.aa_sequence]);
    }

    fn before_delete_hook(&self) -> Result<(), QueryError> {return Ok(());}
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