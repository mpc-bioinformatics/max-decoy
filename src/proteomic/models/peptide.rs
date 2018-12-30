extern crate postgres;

use std::hash::{Hash, Hasher};

use self::postgres::Connection;

use proteomic::models::mass::neutral_loss::NeutralLoss;
use proteomic::models::amino_acids::amino_acid::AminoAcid;
use proteomic::models::persistable::Persistable;

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
            weight: Peptide::calculate_weight(&generalized_aa_sequence),
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

    fn calculate_weight(aa_sequence: &String) -> i64 {
        let mut weight: i64 = NeutralLoss::get("H2O").get_mono_mass();
        for amino_acid_one_letter_code in aa_sequence.as_str().chars() {
            weight += AminoAcid::get(amino_acid_one_letter_code).get_mono_mass();
        }
        return weight;
    }

    pub fn get_weight(&self) -> f64 {
        return self.weight as f64 / WEIGHT_CONVERT_FACTOR;
    }

    pub fn get_aa_sequence(&self) -> &String {
        return &self.aa_sequence;
    }
}

impl Persistable<Peptide, i64, String> for Peptide {
    fn from_sql_row(row: &postgres::rows::Row) -> Result<Self, String> {
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

    fn get_primary_key(&self) -> i64 {
        return self.id;
    }

    fn find(conn: &Connection, primary_key: &i64) -> Result<Self, String> {
        match conn.query("SELECT * FROM peptides WHERE id = $1 LIMIT 1", &[primary_key]) {
            Ok(ref rows) if rows.len() > 0 => Self::from_sql_row(&rows.get(0)),
            Ok(_rows) => Err("NOHIT".to_owned()),
            Err(err) => Err(err.code().unwrap().code().to_owned())
        }
    }

    fn find_by_unique_identifier(conn: &Connection, unique_identifier: &String) -> Result<Self, String> {
        let generalized_aa_sequence: String = Self::gerneralize_aa_sequence(unique_identifier);
        match conn.query(
            "SELECT * FROM peptides WHERE aa_sequence = $1 LIMIT 1",
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
                        return Ok(());
                    },
                    Err(err) => Err(format!("cannot insert nor find peptide '{}'\n\torigina error: {}", self.aa_sequence, err))
                }
            }
            Err(err) => Err(err.code().unwrap().code().to_owned())
        }
    }

    fn update(&mut self, conn: &postgres::Connection) -> Result<(), String> {
        match conn.query(
            Self::get_update_query(),
            &[&self.id, &self.aa_sequence, &self.digest_enzym, &self.number_of_missed_cleavages, &self.weight, &self.length]
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
            return Err("Peptide is not persisted".to_owned());
        }
        match conn.execute("DELETE FROM peptides WHERE id = $1;", &[&self.id]) {
            Ok(_) => {
                self.id = 0;
                return Ok(());
            },
            Err(err) => Err(format!("could not delete Peptide from database; postgresql error is: {}", err))
        }
    }

    fn delete_all(conn: &postgres::Connection) -> Result<(), String> {
        match conn.execute("DELETE FROM peptides WHERE id IS NOT NULL;", &[]) {
            Ok(_) => Ok(()),
            Err(err) => Err(format!("could not delete Peptides from database; postgresql error is: {}", err))
        }
    }


    fn get_table_name() -> &'static str {
        return "peptides";
    }

    fn get_select_primary_key_by_unique_identifier_query() -> &'static str {
        return "SELECT id FROM peptides WHERE aa_sequence = $1 LIMIT 1";
    }

    fn get_insert_query() -> &'static str {
        return "INSERT INTO peptides (aa_sequence, digest_enzym, number_of_missed_cleavages, weight, length) VALUES ($1, $2, $3, $4, $5) ON CONFLICT (aa_sequence, weight) DO NOTHING RETURNING id";
    }

    fn get_update_query() -> &'static str {
        return "UPDATE peptides SET aa_sequence = $2, digest_enzym = $3, number_of_missed_cleavages = $4, weight = $5, length = $6 WHERE id = $1";
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
        match prepared_statement.query(&[&self.aa_sequence, &self.digest_enzym, &self.number_of_missed_cleavages, &self.weight, &self.length]) {
            Ok(ref rows) if rows.len() > 0 => {
                self.id = rows.get(0).get(0);
                return Ok(());
            },
            Ok(_rows) => Err("NORET".to_owned()),
            Err(err) => Err(err.code().unwrap().code().to_owned())
        }
    }

    fn exec_update_statement(&mut self, prepared_statement: &postgres::stmt::Statement) -> Result<(), String> {
        match prepared_statement.query(&[&self.id, &self.aa_sequence, &self.digest_enzym, &self.number_of_missed_cleavages, &self.weight, &self.length]) {
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
