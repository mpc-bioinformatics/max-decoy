extern crate postgres;

use proteomic::models::persistable::{Persistable, QueryError};
use proteomic::models::decoys::base_decoy::BaseDecoy;
use proteomic::models::decoys::modified_decoy::ModifiedDecoy;

pub trait Decoy {
    fn to_string(&self) -> String;
    fn get_header(&self) -> String;
    fn get_aa_sequence(&self) -> String;
    fn get_weight(&self) -> i64;
    fn get_length(&self) -> i32;
    // returns the last char of aa_sequence
    // or '_' if aa_sequence is empty
    fn get_c_terminus_amino_acid(&self) -> char;
    // returns the first char of aa_sequence
    // or '_' if aa_sequence is empty
    fn get_n_terminus_amino_acid(&self) -> char;
    // returns the n-th char of aa_sequence
    // or '_' if idx is larger the aa_sequence
    fn get_amino_acid_at(&self, idx: usize) -> char;

    fn as_fasta_entry(&self) -> String {
        return format!("{}\n{}", self.get_header(), self.get_aa_sequence());
    }
}


pub struct PlainDecoy {
    header: String,
    aa_sequence: String,
    weight: i64
}

impl PlainDecoy {
    pub fn new(header: &str, aa_sequence: &str, weight: &i64) -> Self {
        return Self {
            header: header.to_owned(),
            aa_sequence: aa_sequence.to_owned(),
            weight: *weight
        }
    }

    pub fn count_where_mass_tolerance(conn: &postgres::Connection, lower_mass_limit: i64, upper_mass_limit: i64) -> Result<i64, QueryError>  {
        let condition: &str = "weight BETWEEN $1 AND $2";
        let mut number_of_decoys: i64 = match BaseDecoy::count_where(conn, condition, &[&lower_mass_limit, &upper_mass_limit]) {
            Ok(count) => count,
            Err(err) => return Err(err)
        };
        number_of_decoys += match ModifiedDecoy::count_where(conn, condition, &[&lower_mass_limit, &upper_mass_limit]) {
            Ok(count) => count,
            Err(err) => return Err(err)
        };
        return Ok(number_of_decoys);
    }

    pub fn where_mass_tolerance(conn: &postgres::Connection, lower_mass_limit: i64, upper_mass_limit: i64) -> Result<Vec<Self>, QueryError> {
        let mut decoys: Vec<Self> = match BaseDecoy::find_where_as_prepared_decoys(conn, "weight BETWEEN $1 AND $2", &[&lower_mass_limit, &upper_mass_limit]) {
            Ok(decoys) => decoys,
            Err(err) => return Err(err)
        };
        let condition: String = format!("{}.weight BETWEEN $1 AND $2", ModifiedDecoy::get_table_name());
        let mut decoys_from_modified_decoys: Vec<Self> = match ModifiedDecoy::find_where_as_prepared_decoys(conn, condition.as_str(), &[&lower_mass_limit, &upper_mass_limit]) {
            Ok(decoys) => decoys,
            Err(err) => return Err(err)
        };
        decoys.append(&mut decoys_from_modified_decoys);
        return Ok(decoys);
    }
}

impl Decoy for PlainDecoy {
    fn to_string(&self) -> String {
        return format!(
            "proteomic::models::decoys::decoy::PlainDecoy\n\theader => {}\n\taa_sequence => {}\n\tweight => {}",
            self.header,
            self.aa_sequence,
            self.weight
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