extern crate postgres;

use std::hash::{Hash, Hasher};

use self::postgres::Connection;

use proteomic::models::amino_acids::amino_acid::AminoAcid;
use proteomic::models::persistable::{handle_postgres_error, Persistable, QueryError, QueryOk, FromSqlRowError};

const WEIGHT_CONVERT_FACTOR: f64 = 1000000.0;

/*
 * attributes id, length, number_of_missed_cleavages and weight should be unsigned, but postgresql crate and database does not support it
 * comments behind attributes show databse type
 */
pub struct Peptide {
    id: i64,                            // BIGSERIAL
    aa_sequence: String,                // TEXT
    digest_enzym: String,               // CHAR(5)
    length: i32,                        // INTEGER
    number_of_missed_cleavages: i16,    // SMALLINT
    weight: i64                         // BIGINT
}

impl Peptide {
    pub fn new(aa_sequence: String, digest_enzym: String, number_of_missed_cleavages: i16) -> Peptide {
        let generalized_aa_sequence: String = Peptide::gerneralize_aa_sequence(&aa_sequence);
        return Peptide{
            id: 0,
            length: generalized_aa_sequence.len() as i32,
            weight: AminoAcid::get_sequence_weight(&generalized_aa_sequence),
            aa_sequence: generalized_aa_sequence,
            digest_enzym: digest_enzym,
            number_of_missed_cleavages: number_of_missed_cleavages
        }
    }

    pub fn is_new(&self) -> bool {
        return self.id < 1;
    }

    pub fn to_string(&self) -> String {
        return format!("{}: {}\n\tdigested with => {}\n\tmissed_cleavages => {}\n\tweight => {}\n\tlength => {}", self.id, self.aa_sequence, self.digest_enzym, self.number_of_missed_cleavages, self.weight, self.length);
    }

    pub fn gerneralize_aa_sequence(aa_sequence: &String) -> String {
        return aa_sequence.replace("I", "J").replace("L", "J");
    }

    pub fn get_weight(&self) -> f64 {
        return self.weight as f64 / WEIGHT_CONVERT_FACTOR;
    }

    pub fn get_aa_sequence(&self) -> &String {
        return &self.aa_sequence;
    }
}

impl Persistable<Peptide, i64, String> for Peptide {
    fn from_sql_row(row: &postgres::rows::Row) -> Result<Self, FromSqlRowError> {
        return Ok(
            Self{
                id: row.get(0),
                aa_sequence: row.get(1),
                length: row.get(2),
                number_of_missed_cleavages: row.get(3),
                weight: row.get(4),
                digest_enzym: row.get(5)
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
        return "peptides";
    }

    fn is_persisted(&self) -> bool {
        return self.id > 0;
    }

    fn find_query() -> &'static str {
        return "SELECT * FROM peptides WHERE id = $1 LIMIT 1;";
    }

    fn create_query() -> &'static str {
        return "INSERT INTO peptides (aa_sequence, digest_enzym, number_of_missed_cleavages, weight, length) VALUES ($1, $2, $3, $4, $5) ON CONFLICT (aa_sequence, weight) DO NOTHING RETURNING id;";
    }

    fn create_attributes(&self) -> Box<Vec<&postgres::types::ToSql>>{
        return Box::new(vec![&self.aa_sequence, &self.digest_enzym, &self.number_of_missed_cleavages, &self.weight, &self.length]);
    }

    fn update_query() -> &'static str{
        return "UPDATE peptides SET aa_sequence = $2, digest_enzym = $3, number_of_missed_cleavages = $4, weight = $5, length = $6 WHERE id = $1;";
    }

    fn update_attributes(&self) -> Box<Vec<&postgres::types::ToSql>>{
        return Box::new(vec![&self.id, &self.aa_sequence, &self.digest_enzym, &self.number_of_missed_cleavages, &self.weight, &self.length]);
    }

    fn delete_query() -> &'static str {
        return "DELETE FROM peptides WHERE id = $1;";
    }

    fn delete_attributes(&self) -> Box<Vec<&postgres::types::ToSql>> {
        return Box::new(vec![&self.id]);
    }

    fn delete_all_query() -> &'static str {
        return "DELETE FROM peptides WHERE id IS NOT NULL;";
    }

    fn exists_query() -> &'static str {
        return "SELECT id FROM peptides WHERE aa_sequence = $1 LIMIT 1;";
    }

    fn exists_attributes(&self) -> Box<Vec<&postgres::types::ToSql>> {
        return Box::new(vec![&self.aa_sequence]);
    }

    fn before_delete_hook(&self) -> Result<(), QueryError> {return Ok(());}
}

// PartialEq-implementation to use this type in a HashSet
impl PartialEq for Peptide {
    fn eq(&self, other: &Peptide) -> bool {
       return self.aa_sequence == *other.get_aa_sequence();
    }
}

// Eq-implementation to use this type in a HashSet
impl Eq for Peptide {}

// Hash-implementation to use this type in a HashSet
impl Hash for Peptide {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.aa_sequence.hash(state);
    }
}

impl Clone for Peptide {
    fn clone(&self) -> Peptide {
        return Peptide{
            id: self.id,
            aa_sequence:  String::from(self.aa_sequence.as_str()),
            length: self.length,
            number_of_missed_cleavages: self.number_of_missed_cleavages,
            weight: self.weight,
            digest_enzym: String::from(self.digest_enzym.as_str()),
        }
    }
}
