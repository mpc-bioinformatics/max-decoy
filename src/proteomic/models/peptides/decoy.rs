use std::hash::{Hash, Hasher};
use std::collections::HashMap;

use proteomic::models::mass;
use proteomic::models::persistable::{Persistable, QueryOk, QueryError, FromSqlRowError};
use proteomic::models::peptides::peptide_interface::PeptideInterface;
use proteomic::models::peptides::peptide::Peptide;
use proteomic::models::amino_acids::amino_acid::AminoAcid;

pub const DECOY_HEADER_START: &'static str = ">DECOY_";

/// In fact, this is the nearly same implementation as Peptide.
/// Major difference is the table name. There should be implemented a trait for this Peptide/Decoy at a given time.
pub struct Decoy {
    id: i64,                                // BIGSERIAL
    aa_sequence: String,                    // TEXT
    length: i32,                            // INTEGER
    number_of_missed_cleavages: i16,        // SMALLINT
    weight: i64,                            // BIGINT
    amino_acids_counts: HashMap<char, i16>,  // columns <one_letter_code>_count, type SMALLINT
    modification_summary: String
}

impl Decoy {
    pub fn new(aa_sequence: &str, number_of_missed_cleavages: i16) -> Self {
        let generalized_aa_sequence: String = AminoAcid::gerneralize_sequence(aa_sequence);
        return Self {
            id: 0,
            length: generalized_aa_sequence.len() as i32,
            weight: AminoAcid::get_sequence_weight(generalized_aa_sequence.as_str()),
            amino_acids_counts: *Self::count_amino_acids(generalized_aa_sequence.as_str()),
            aa_sequence: generalized_aa_sequence,
            number_of_missed_cleavages: number_of_missed_cleavages,
            modification_summary: String::new()
        }
    }

    pub fn get_count_for_amino_acid(&self, amino_acid_one_letter_code: &char) -> &i16 {
        match self.amino_acids_counts.get(amino_acid_one_letter_code) {
            Some(count) => count,
            None => &0,
        }
    }

    pub fn get_header(&self) -> String {
        let mut header = format!("{}{}", DECOY_HEADER_START, self.aa_sequence.as_str());
        if self.modification_summary.len() > 0 {
            header.push_str(format!(" ModRes={}", self.modification_summary).as_str());
        }
        return header;
    }

    pub fn set_modification_summary(&mut self, modification_summary: &str) {
        self.modification_summary = modification_summary.to_owned();
    }

    pub fn is_peptide(&self, conn: &postgres::Connection) -> bool {
        match Peptide::exists_where(conn, "aa_sequence = $1", &[&self.aa_sequence]) {
            Ok(query_ok) => match query_ok {
                QueryOk::Exists => return true,
                _ => panic!("proteomic::models::peptides::decoy::Decoy.is_peptide(): Got QueryOk which should not be returned")
            },
            Err(query_err) => match query_err {
                QueryError::NoMatch => return false,
                _ => panic!("proteomic::models::peptides::decoy::Decoy.is_peptide(): Err: {}", query_err)
            }
        }
    }
}

impl PeptideInterface for Decoy {
    fn to_string(&self) -> String {
        return format!(
            "proteomic::models::peptides::decoy::Decoy\n\tid => {}\n\taa_sequence => {}\n\tlength => {}\n\tnumber_of_missed_cleavages => {}\n\tweight => {}\n\taa_counts => {:?}",
            self.id,
            self.aa_sequence,
            self.length,
            self.number_of_missed_cleavages,
            mass::convert_mass_to_float(self.weight),
            self.amino_acids_counts
        );
    }

    fn get_aa_sequence(&self) -> &str {
        return self.aa_sequence.as_str();
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

impl Persistable<Decoy, i64, String> for Decoy {
    fn from_sql_row(row: &postgres::rows::Row) -> Result<Self, FromSqlRowError> {
        return Ok(
            Self{
                id: row.get(0),
                aa_sequence: row.get::<usize, String>(1).trim().to_owned(),     // remember that aa_sequence is saved as CHAR(60) which means the database will fill it up with whitespaces if it is shorter
                length: row.get(2),
                number_of_missed_cleavages: row.get(3),
                weight: row.get(4),
                amino_acids_counts: [
                    ('r', row.get(5)),
                    ('n', row.get(6)),
                    ('d', row.get(7)),
                    ('c', row.get(8)),
                    ('e', row.get(9)),
                    ('q', row.get(10)),
                    ('g', row.get(11)),
                    ('h', row.get(12)),
                    ('j', row.get(13)),
                    ('k', row.get(14)),
                    ('m', row.get(15)),
                    ('f', row.get(16)),
                    ('p', row.get(17)),
                    ('o', row.get(18)),
                    ('s', row.get(19)),
                    ('t', row.get(20)),
                    ('u', row.get(21)),
                    ('v', row.get(22)),
                    ('w', row.get(23)),
                    ('y', row.get(24))
                ].iter().cloned().collect()
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
        return "decoys";
    }

    fn is_persisted(&self) -> bool {
        return self.id > 0;
    }

    fn find_query() -> &'static str {
        return "SELECT * FROM decoys WHERE id = $1 LIMIT 1;";
    }

    fn create_query() -> &'static str {
        return "INSERT INTO decoys (aa_sequence, number_of_missed_cleavages, weight, length, r_count, n_count, d_count, c_count, e_count, q_count, g_count, h_count, j_count, k_count, m_count, f_count, p_count, o_count, s_count, t_count, u_count, v_count, w_count, y_count) VALUES ($1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22, $23, $24) ON CONFLICT (aa_sequence, weight) DO NOTHING RETURNING id;";
    }

    fn create_attributes(&self) -> Box<Vec<&postgres::types::ToSql>>{
        return Box::new(vec![
            &self.aa_sequence,
            &self.number_of_missed_cleavages,
            &self.weight,
            &self.length,
            self.get_count_for_amino_acid(&'R'),
            self.get_count_for_amino_acid(&'N'),
            self.get_count_for_amino_acid(&'D'),
            self.get_count_for_amino_acid(&'C'),
            self.get_count_for_amino_acid(&'E'),
            self.get_count_for_amino_acid(&'Q'),
            self.get_count_for_amino_acid(&'G'),
            self.get_count_for_amino_acid(&'H'),
            self.get_count_for_amino_acid(&'J'),
            self.get_count_for_amino_acid(&'K'),
            self.get_count_for_amino_acid(&'M'),
            self.get_count_for_amino_acid(&'F'),
            self.get_count_for_amino_acid(&'P'),
            self.get_count_for_amino_acid(&'O'),
            self.get_count_for_amino_acid(&'S'),
            self.get_count_for_amino_acid(&'T'),
            self.get_count_for_amino_acid(&'U'),
            self.get_count_for_amino_acid(&'V'),
            self.get_count_for_amino_acid(&'W'),
            self.get_count_for_amino_acid(&'Y')
        ]);
    }

    fn update_query() -> &'static str{
        return "UPDATE decoys SET aa_sequence = $2, number_of_missed_cleavages = $4, weight = $5, length = $6, r_count = $7, n_count = $8, d_count = $9, c_count = $10, e_count = $11, q_count = $12, g_count = $13, h_count = $14, j_count = $15, k_count = $16, m_count = $17, f_count = $18, p_count = $19, o_count = $20, s_count = $21, t_count = $22, u_count = $23, v_count = $24, w_count = $25, y_count = $26 WHERE id = $1;";
    }

    fn update_attributes(&self) -> Box<Vec<&postgres::types::ToSql>>{
        return Box::new(vec![
            &self.id,
            &self.aa_sequence,
            &self.number_of_missed_cleavages,
            &self.weight,
            &self.length,
            self.get_count_for_amino_acid(&'R'),
            self.get_count_for_amino_acid(&'N'),
            self.get_count_for_amino_acid(&'D'),
            self.get_count_for_amino_acid(&'C'),
            self.get_count_for_amino_acid(&'E'),
            self.get_count_for_amino_acid(&'Q'),
            self.get_count_for_amino_acid(&'G'),
            self.get_count_for_amino_acid(&'H'),
            self.get_count_for_amino_acid(&'J'),
            self.get_count_for_amino_acid(&'K'),
            self.get_count_for_amino_acid(&'M'),
            self.get_count_for_amino_acid(&'F'),
            self.get_count_for_amino_acid(&'P'),
            self.get_count_for_amino_acid(&'O'),
            self.get_count_for_amino_acid(&'S'),
            self.get_count_for_amino_acid(&'T'),
            self.get_count_for_amino_acid(&'U'),
            self.get_count_for_amino_acid(&'V'),
            self.get_count_for_amino_acid(&'W'),
            self.get_count_for_amino_acid(&'Y')
        ]);
    }

    fn delete_query() -> &'static str {
        return "DELETE FROM decoys WHERE id = $1;";
    }

    fn delete_attributes(&self) -> Box<Vec<&postgres::types::ToSql>> {
        return Box::new(vec![&self.id]);
    }

    fn delete_all_query() -> &'static str {
        return "DELETE FROM decoys WHERE id IS NOT NULL;";
    }

    fn exists_query() -> &'static str {
        return "SELECT id FROM decoys WHERE aa_sequence = $1 LIMIT 1;";
    }

    fn exists_attributes(&self) -> Box<Vec<&postgres::types::ToSql>> {
        return Box::new(vec![&self.aa_sequence]);
    }

    fn before_delete_hook(&self) -> Result<(), QueryError> {return Ok(());}
}

// PartialEq-implementation to use this type in a HashSet
impl PartialEq for Decoy {
    fn eq(&self, other: &Decoy) -> bool {
        return self.aa_sequence.eq(&other.get_aa_sequence());
    }
}

// Eq-implementation to use this type in a HashSet
impl Eq for Decoy {}

// Hash-implementation to use this type in a HashSet
impl Hash for Decoy {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.aa_sequence.hash(state);
    }
}