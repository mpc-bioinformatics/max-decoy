extern crate postgres;

use std::hash::{Hash, Hasher};
use std::error::Error;

use self::postgres::Connection;

use proteomic::utility::{amino_acid, mass};
use proteomic::models::persistable::Persistable;
use proteomic::models::collection::Collectable;

const WEIGHT_CONVERT_FACTOR: f32 = 1000000.0;

pub struct Peptide {
    id: i32,
    aa_sequence: String,
    digest_enzym: String,
    length: i32,
    number_of_missed_cleavages: i32,
    weight: i32,
    is_persisted: bool
}

impl Peptide {
    pub fn new(aa_sequence: String, digest_enzym: String, number_of_missed_cleavages: i32) -> Peptide {
        let generalized_aa_sequence: String = Peptide::gerneralize_aa_sequence(&aa_sequence);
        return Peptide{
            id: -1,
            length: generalized_aa_sequence.len() as i32,
            weight: Peptide::calculate_weight(&generalized_aa_sequence),
            aa_sequence: generalized_aa_sequence,
            digest_enzym: digest_enzym,
            number_of_missed_cleavages: number_of_missed_cleavages,
            is_persisted: false
        }
    }

    pub fn is_new(&self) -> bool {
        return self.id < 0;
    }

    pub fn to_string(&self) -> String {
        return format!("{}: {}\n\tdigested with => {}\n\tmissed_cleavages => {}\n\tweight => {}\n\tlength => {}", self.id, self.aa_sequence, self.digest_enzym, self.number_of_missed_cleavages, self.weight, self.length);
    }

    pub fn gerneralize_aa_sequence(aa_sequence: &String) -> String {
        return aa_sequence.replace("I", "J").replace("L", "J");
    }

    fn calculate_weight(aa_sequence: &String) -> i32 {
        let mut weight: f32 = 0.0;
        for amino_acid_one_letter_code in aa_sequence.as_str().chars() {
            weight += amino_acid::get(amino_acid_one_letter_code).get_mono_mass();
        }
        weight += mass::get_neutral_loss("H2O").get_mono_mass();
        return (weight * WEIGHT_CONVERT_FACTOR) as i32;
    }

    pub fn get_weight(&self) -> f32 {
        return self.weight as f32 / WEIGHT_CONVERT_FACTOR;
    }

    pub fn get_aa_sequence(&self) -> &String {
        return &self.aa_sequence;
    }
}

impl Persistable<Peptide, i32, String> for Peptide {
    fn get_primary_key(&self) -> i32 {
        return self.id;
    }

    fn find(conn: &Connection, primary_key: &i32) -> Result<Self, String> {
        match conn.query("SELECT * FROM peptides WHERE id = $1 LIMIT 1", &[primary_key]) {
            Ok(rows) =>{
                if rows.len() > 0 {
                    Ok(
                        Peptide{
                            id: rows.get(0).get(0),
                            aa_sequence: rows.get(0).get(1),
                            length: rows.get(0).get(2),
                            number_of_missed_cleavages: rows.get(0).get(3),
                            weight: rows.get(0).get(4),
                            digest_enzym: rows.get(0).get(5),
                            is_persisted: true
                        }
                    )
                } else {
                    Err("peptide not found".to_owned())
                }
            },
            Err(_err) => Err("some error occured: Peptide::find".to_owned())
        }
    }

    fn find_by_unique_identifier(conn: &Connection, unique_identifier: &String) -> Result<Self, String> {
        let generalized_aa_sequence: String = Self::gerneralize_aa_sequence(unique_identifier);
        match conn.query(
            "SELECT * FROM peptides WHERE aa_sequence = $1 LIMIT 1",
            &[&generalized_aa_sequence]
        ) {
            Ok(ref rows) if rows.len() > 0 => Ok(
                Peptide{
                    id: rows.get(0).get(0),
                    aa_sequence: rows.get(0).get(1),
                    length: rows.get(0).get(2),
                    number_of_missed_cleavages: rows.get(0).get(3),
                    weight: rows.get(0).get(4),
                    digest_enzym: rows.get(0).get(5),
                    is_persisted: true
                }
            ),
            Ok(_rows) => Err("peptide not found".to_owned()),
            Err(err) => Err(err.description().to_owned())
        }
    }


    fn create(&mut self, conn: &postgres::Connection) -> Result<(), String> {
        match conn.query(
            Self::get_insert_query(),
            &[&self.aa_sequence, &self.digest_enzym, &self.number_of_missed_cleavages, &self.weight, &self.length]
        ) {
            Ok(ref rows) if rows.len() > 0 => {
                self.id =  rows.get(0).get(0);
                return Ok(());
            },
            Ok(_rows) => {
                // zero rows means there are a conflict on update, so the peptides exists already
                match Self::find_by_unique_identifier(conn, &self.aa_sequence) {
                    Ok(peptide) => {
                        self.id = peptide.get_primary_key();
                        self.is_persisted = true;
                        return Ok(());
                    },
                    Err(err) => Err(format!("cannot insert nor find peptide '{}'\n\torigina error: {}", self.aa_sequence, err))
                }
            }
            Err(err) => Err(err.description().to_owned())
        }
    }

    fn update(&mut self, conn: &postgres::Connection) -> Result<(), String> {
        match conn.query(
            Self::get_update_query(),
            &[&self.id, &self.aa_sequence, &self.digest_enzym, &self.number_of_missed_cleavages, &self.weight, &self.length]
        ) {
            Ok(ref rows) if rows.len() > 0 => Ok(()),
            Ok(_rows) => Err("updateing peptide does not return anything".to_owned()),
            Err(err) => Err(err.description().to_owned())
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
        return "SELECT id FROM peptides WHERE aa_sequence = $1 LIMIT 1";
    }

    fn get_insert_query() -> &'static str {
        return "INSERT INTO peptides (aa_sequence, digest_enzym, number_of_missed_cleavages, weight, length) VALUES ($1, $2, $3, $4, $5) ON CONFLICT (aa_sequence) DO NOTHING RETURNING id";
    }

    fn get_update_query() -> &'static str {
        return "UPDATE peptides SET aa_sequence = $2, digest_enzym = $3, number_of_missed_cleavages = $4, weight = $5, length = $6 WHERE id = $1";
    }

    fn exec_select_primary_key_by_unique_identifier_statement(&mut self, prepared_statement: &postgres::stmt::Statement) -> Result<(), String> {
        match prepared_statement.query(&[&self.aa_sequence]) {
            Ok(ref rows) if rows.len() > 0 => {
                self.id = rows.get(0).get(0);
                self.is_persisted = true;
                return Ok(());
            },
            Ok(_rows) => Err("peptides not found".to_owned()),
            Err(err) => Err(err.description().to_owned())
        }
    }

    fn exec_insert_statement(&mut self, prepared_statement: &postgres::stmt::Statement) -> Result<(), String> {
        match prepared_statement.query(&[&self.aa_sequence, &self.digest_enzym, &self.number_of_missed_cleavages, &self.weight, &self.length]) {
            Ok(ref rows) if rows.len() > 0 => {
                self.id = rows.get(0).get(0);
                self.is_persisted = true;
                return Ok(());
            },
            Ok(_rows) => Err("inserting peptides does not return anything".to_owned()),
            Err(err) => Err(err.description().to_owned())
        }
    }

    fn exec_update_statement(&mut self, prepared_statement: &postgres::stmt::Statement) -> Result<(), String> {
        match prepared_statement.query(&[&self.id, &self.aa_sequence, &self.digest_enzym, &self.number_of_missed_cleavages, &self.weight, &self.length]) {
            Ok(ref rows) if rows.len() > 0 => Ok(()),
            Ok(_rows) => Err("updating peptides does not return anything".to_owned()),
            Err(err) => Err(err.description().to_owned())
        }
    }


    fn is_persisted(&self) -> bool {
        return self.is_persisted;
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

impl Clone for Peptide {
    fn clone(&self) -> Peptide {
        return Peptide{
            id: self.id,
            aa_sequence:  String::from(self.aa_sequence.as_str()),
            length: self.length,
            number_of_missed_cleavages: self.number_of_missed_cleavages,
            weight: self.weight,
            digest_enzym: String::from(self.digest_enzym.as_str()),
            is_persisted: self.is_persisted
        }
    }
}