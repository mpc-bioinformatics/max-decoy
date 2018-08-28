use std::hash::{Hash, Hasher};

use super::postgres::Connection;
use super::postgres::rows::Rows;
use super::postgres::Result;

use proteomic::utility::amino_acid;
use proteomic::models::persistable::Persistable;
use proteomic::models::collection::Collectable;
use super::postgres::stmt::Statement;

pub struct Peptide {
    id: i32,
    aa_sequence: String,
    digest_enzym: String,
    length: i32,
    position: usize,
    number_of_missed_cleavages: i32,
    weight: f64
}

impl Peptide {
    pub fn new(aa_sequence: String, digest_enzym: String, position: usize, number_of_missed_cleavages: i32) -> Peptide {
        return Peptide{
            id: -1,
            length: aa_sequence.len() as i32,
            weight: Peptide::calculate_weight(&aa_sequence),
            aa_sequence: aa_sequence,
            digest_enzym: digest_enzym,
            position: position,
            number_of_missed_cleavages: number_of_missed_cleavages
        }
    }

    pub fn print(&self) {
        println!("{}\n\tdigested with => {}\n\tpos => {}\n\tmissed_cleavages => {}\n\tweight => {}\n\tlength => {}", self.aa_sequence, self.digest_enzym, self.position, self.number_of_missed_cleavages, self.weight, self.length);
    }

    fn calculate_weight(aa_sequence: &String) -> f64 {
        let mut weight: f64 = 0.0;
        for amino_acid_one_letter_code in aa_sequence.as_str().chars() {
            weight += amino_acid::get(amino_acid_one_letter_code).get_mono_mass();
        }
        return weight;
    }

    pub fn create(&mut self) {
        let conn = super::get_db_connection();
        for row in self.insert_query(&conn).unwrap().iter() {
            let id: i32 = row.get("id");
            self.id = id;
            break;
        }
    }

    pub fn update(&self) {
        let conn = super::get_db_connection();
        self.update_query(&conn);
    }

    pub fn save(&mut self) {
        if self.id > 0 {
            self.update();
        } else {
            self.create();
        }
    }
}

impl Persistable for Peptide {
    fn get_id(&self) -> i32 {
        return self.id;
    }
    
    fn get_insert_statement() -> &'static str {
        return "INSERT INTO peptides (aa_sequence, digest_enzym, number_of_missed_cleavages, weight, length) VALUES ($1, $2, $3, $4, $5) RETURNING id";
    }

    fn insert_query(&self, connection: &Connection) -> Result<Rows> {
        return connection.query(
            Peptide::get_insert_statement(),
            &[&self.aa_sequence, &self.digest_enzym, &self.number_of_missed_cleavages, &self.weight]
        );
    }

    fn insert_statement(&self, prepared_statement: &Statement) {
        prepared_statement.execute(&[&self.aa_sequence, &self.digest_enzym, &self.number_of_missed_cleavages, &self.weight, &self.length]);
    }

    fn get_update_statement() -> &'static str {
        return "UPDATE peptides SET aa_sequence = $2, digest_enzym = $3, number_of_missed_cleavages = $4, weight = $5, length = $6 WHERE id = $1";
    }

    fn update_query(&self, connection: &Connection) -> Result<Rows> {
        return connection.query(
            Peptide::get_update_statement(),
            &[&self.id, &self.aa_sequence, &self.digest_enzym, &self.number_of_missed_cleavages, &self.weight]
        );
    }

    fn update_statement(&self, prepared_statement: &Statement) {
        prepared_statement.execute(&[&self.id, &self.aa_sequence, &self.digest_enzym, &self.number_of_missed_cleavages, &self.weight, &self.length]);
    }
}

impl Collectable for Peptide {
    fn get_collection_identifier(&self) -> &String {
        return &self.aa_sequence;
    }
}

// PartialEq-implementation to use this type in a HashSet
impl PartialEq for Peptide {
    fn eq(&self, other: &Peptide) -> bool {
       return self.aa_sequence.eq(&other.aa_sequence);
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