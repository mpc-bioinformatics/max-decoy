extern crate postgres;

use std::hash::{Hash, Hasher};

use self::postgres::Connection;
use self::postgres::rows::Rows;
use self::postgres::Result;
use self::postgres::stmt::Statement;

use proteomic::utility::amino_acid;
use proteomic::models::persistable::Persistable;
use proteomic::models::collection::Collectable;


pub struct Peptide {
    id: i32,
    aa_sequence: String,
    digest_enzym: String,
    length: i32,
    position: usize,
    number_of_missed_cleavages: i32,
    weight: i32
}

impl Peptide {
    pub fn new(aa_sequence: String, digest_enzym: String, position: usize, number_of_missed_cleavages: i32) -> Peptide {
        let generalized_aa_sequence: String = Peptide::gerneralize_aa_sequence(&aa_sequence);
        return Peptide{
            id: -1,
            length: generalized_aa_sequence.len() as i32,
            weight: Peptide::calculate_weight(&generalized_aa_sequence),
            aa_sequence: generalized_aa_sequence,
            digest_enzym: digest_enzym,
            position: position,
            number_of_missed_cleavages: number_of_missed_cleavages
        }
    }

    pub fn print(&self) {
        println!("{}\n\tdigested with => {}\n\tpos => {}\n\tmissed_cleavages => {}\n\tweight => {}\n\tlength => {}", self.aa_sequence, self.digest_enzym, self.position, self.number_of_missed_cleavages, self.weight, self.length);
    }

    fn gerneralize_aa_sequence(aa_sequence: &String) -> String {
        return aa_sequence.replace("I", "J").replace("L", "J");
    }

    fn calculate_weight(aa_sequence: &String) -> i32 {
        let mut weight: f32 = 0.0;
        for amino_acid_one_letter_code in aa_sequence.as_str().chars() {
            weight += amino_acid::get(amino_acid_one_letter_code).get_mono_mass();
        }
        return (weight * 100000.0) as i32;
    }

    pub fn get_weight(&self) -> f32 {
        return self.weight as f32 / 1000000.0;
    }

    pub fn get_aa_sequence(&self) -> &String {
        return &self.aa_sequence;
    }

    pub fn create(&mut self, conn: &Connection) {
        for row in self.execute_insert_query(conn).unwrap().iter() {
            let id: i32 = row.get("id");
            self.id = id;
            break;
        }
    }

    pub fn update(&self, conn: &Connection) {
        self.execute_update_query(conn);
    }

    pub fn save(&mut self, conn: &Connection) {
        if self.id > 0 {
            self.update(conn);
        } else {
            self.create(conn);
        }
    }

    pub fn exists(&self, conn: &Connection) -> bool {
        for row in self.exists_query(conn).unwrap().iter() {
            return row.get::<usize, bool>(0);
        }
        return false;
    }
}

impl Persistable for Peptide {
    fn get_id(&self) -> i32 {
        return self.id;
    }

    fn get_insert_statement() -> &'static str {
        return "INSERT INTO peptides (aa_sequence, digest_enzym, number_of_missed_cleavages, weight, length) VALUES ($1, $2, $3, $4, $5) ON CONFLICT (aa_sequence) DO NOTHING";
    }

    fn execute_insert_query(&self, connection: &Connection) -> Result<Rows> {
        return connection.query(
            Peptide::get_insert_statement(),
            &[&self.aa_sequence, &self.digest_enzym, &self.number_of_missed_cleavages, &self.weight, &self.length]
        );
    }

    fn execute_insert_statement(&self, prepared_statement: &Statement) {
        prepared_statement.execute(&[&self.aa_sequence, &self.digest_enzym, &self.number_of_missed_cleavages, &self.weight, &self.length]);
    }

    fn get_update_statement() -> &'static str {
        return "UPDATE peptides SET aa_sequence = $2, digest_enzym = $3, number_of_missed_cleavages = $4, weight = $5, length = $6 WHERE id = $1";
    }

    fn execute_update_query(&self, connection: &Connection) -> Result<Rows> {
        return connection.query(
            Peptide::get_update_statement(),
            &[&self.id, &self.aa_sequence, &self.digest_enzym, &self.number_of_missed_cleavages, &self.weight, &self.length]
        );
    }

    fn execute_update_statement(&self, prepared_statement: &Statement) {
        prepared_statement.execute(&[&self.id, &self.aa_sequence, &self.digest_enzym, &self.number_of_missed_cleavages, &self.weight, &self.length]);
    }

    fn get_exists_statement() -> &'static str {
        return "SELECT EXISTS(SELECT 1 FROM peptides WHERE aa_sequence = $1)";
    }

    fn exists_query(&self, connection: &Connection) -> Result<Rows> {
        return connection.query(
            Peptide::get_exists_statement(),
            &[&self.aa_sequence]
        );
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