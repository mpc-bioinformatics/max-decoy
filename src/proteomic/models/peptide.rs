use std::hash::{Hash, Hasher};
use std::collections::HashMap;

use proteomic::models::mass;
use proteomic::models::amino_acids::amino_acid::AminoAcid;
use proteomic::models::persistable::{Persistable, QueryError, FromSqlRowError};

pub const AMINO_ACIDS_FOR_COUNTING: &'static [char] = &['R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'J', 'K', 'M', 'F', 'P', 'O', 'S', 'T', 'U', 'V', 'W', 'Y'];

/*
 * attributes id, length, number_of_missed_cleavages and weight should be unsigned, but postgresql crate and database does not support it
 * comments behind attributes show databse type
 */
pub struct Peptide {
    id: i64,                                // BIGSERIAL
    aa_sequence: String,                    // TEXT
    length: i32,                            // INTEGER
    number_of_missed_cleavages: i16,        // SMALLINT
    weight: i64,                            // BIGINT
    amino_acids_counts: HashMap<char, i16>  // columns <one_letter_code>_count, type SMALLINT
}

impl Peptide {
    pub fn new(aa_sequence: &str, number_of_missed_cleavages: i16) -> Peptide {
        let generalized_aa_sequence: String = AminoAcid::gerneralize_sequence(aa_sequence);
        let mut amino_acids_counts: HashMap<char, i16> = HashMap::new();
        for one_letter_code in AMINO_ACIDS_FOR_COUNTING {
            amino_acids_counts.insert(*one_letter_code, generalized_aa_sequence.matches(*one_letter_code).count() as i16);
        }
        return Peptide{
            id: 0,
            length: generalized_aa_sequence.len() as i32,
            weight: AminoAcid::get_sequence_weight(&generalized_aa_sequence),
            aa_sequence: generalized_aa_sequence,
            number_of_missed_cleavages: number_of_missed_cleavages,
            amino_acids_counts: amino_acids_counts
        }
    }

    pub fn to_string(&self) -> String {
        return format!(
            "proteomic::models::peptide::Peptide\n\tid => {}\n\taa_sequence => {}\n\tlength => {}\n\tnumber_of_missed_cleavages => {}\n\tweight => {}",
            self.id,
            self.aa_sequence,
            self.length,
            self.number_of_missed_cleavages,
            mass::convert_mass_to_float(self.weight)
        );
    }

    pub fn get_weight(&self) -> i64 {
        return self.weight;
    }

    pub fn get_aa_sequence(&self) -> &String {
        return &self.aa_sequence;
    }

    pub fn get_count_for_amino_acid(&self, amino_acid_one_letter_code: &char) -> &i16 {
        match self.amino_acids_counts.get(amino_acid_one_letter_code) {
            Some(count) => count,
            None => &0,
        }
    }
}

impl Persistable<Peptide, i64, String> for Peptide {
    fn from_sql_row(row: &postgres::rows::Row) -> Result<Self, FromSqlRowError> {
        return Ok(
            Self{
                id: row.get(0),
                aa_sequence: row.get::<usize, String>(1).trim().to_owned(),     // remember that aa_sequence is saved as CHAR(60) which means the database will fill it up with whitespaces if it is shorter
                length: row.get(2),
                number_of_missed_cleavages: row.get(3),
                weight: row.get(4),
                amino_acids_counts: [
                    ('r', row.get(6)),
                    ('n', row.get(7)),
                    ('d', row.get(8)),
                    ('c', row.get(9)),
                    ('e', row.get(10)),
                    ('q', row.get(11)),
                    ('g', row.get(12)),
                    ('h', row.get(13)),
                    ('j', row.get(14)),
                    ('k', row.get(15)),
                    ('m', row.get(16)),
                    ('f', row.get(17)),
                    ('p', row.get(18)),
                    ('o', row.get(19)),
                    ('s', row.get(20)),
                    ('t', row.get(21)),
                    ('u', row.get(22)),
                    ('v', row.get(23)),
                    ('w', row.get(24)),
                    ('y', row.get(25))
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
        return "peptides";
    }

    fn is_persisted(&self) -> bool {
        return self.id > 0;
    }

    fn find_query() -> &'static str {
        return "SELECT * FROM peptides WHERE id = $1 LIMIT 1;";
    }

    fn create_query() -> &'static str {
        return "INSERT INTO peptides (aa_sequence, number_of_missed_cleavages, weight, length, r_count, n_count, d_count, c_count, e_count, q_count, g_count, h_count, j_count, k_count, m_count, f_count, p_count, o_count, s_count, t_count, u_count, v_count, w_count, y_count) VALUES ($1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22, $23, $24, $25) ON CONFLICT (aa_sequence, weight) DO NOTHING RETURNING id;";
    }

    fn create_attributes(&self) -> Box<Vec<&postgres::types::ToSql>>{
        return Box::new(vec![
            &self.aa_sequence,
            &self.number_of_missed_cleavages,
            &self.weight,
            &self.length,
            self.get_count_for_amino_acid(&'r'),
            self.get_count_for_amino_acid(&'n'),
            self.get_count_for_amino_acid(&'d'),
            self.get_count_for_amino_acid(&'c'),
            self.get_count_for_amino_acid(&'e'),
            self.get_count_for_amino_acid(&'q'),
            self.get_count_for_amino_acid(&'g'),
            self.get_count_for_amino_acid(&'h'),
            self.get_count_for_amino_acid(&'j'),
            self.get_count_for_amino_acid(&'k'),
            self.get_count_for_amino_acid(&'m'),
            self.get_count_for_amino_acid(&'f'),
            self.get_count_for_amino_acid(&'p'),
            self.get_count_for_amino_acid(&'o'),
            self.get_count_for_amino_acid(&'s'),
            self.get_count_for_amino_acid(&'t'),
            self.get_count_for_amino_acid(&'u'),
            self.get_count_for_amino_acid(&'v'),
            self.get_count_for_amino_acid(&'w'),
            self.get_count_for_amino_acid(&'y')
        ]);
    }

    fn update_query() -> &'static str{
        return "UPDATE peptides SET aa_sequence = $2, number_of_missed_cleavages = $4, weight = $5, length = $6, r_count = $7, n_count = $8, d_count = $9, c_count = $10, e_count = $11, q_count = $12, g_count = $13, h_count = $14, j_count = $15, k_count = $16, m_count = $17, f_count = $18, p_count = $19, o_count = $20, s_count = $21, t_count = $22, u_count = $23, v_count = $24, w_count = $25, y_count = $26 WHERE id = $1;";
    }

    fn update_attributes(&self) -> Box<Vec<&postgres::types::ToSql>>{
        return Box::new(vec![
            &self.id,
            &self.aa_sequence,
            &self.number_of_missed_cleavages,
            &self.weight,
            &self.length,
            self.get_count_for_amino_acid(&'r'),
            self.get_count_for_amino_acid(&'n'),
            self.get_count_for_amino_acid(&'d'),
            self.get_count_for_amino_acid(&'c'),
            self.get_count_for_amino_acid(&'e'),
            self.get_count_for_amino_acid(&'q'),
            self.get_count_for_amino_acid(&'g'),
            self.get_count_for_amino_acid(&'h'),
            self.get_count_for_amino_acid(&'j'),
            self.get_count_for_amino_acid(&'k'),
            self.get_count_for_amino_acid(&'m'),
            self.get_count_for_amino_acid(&'f'),
            self.get_count_for_amino_acid(&'p'),
            self.get_count_for_amino_acid(&'o'),
            self.get_count_for_amino_acid(&'s'),
            self.get_count_for_amino_acid(&'t'),
            self.get_count_for_amino_acid(&'u'),
            self.get_count_for_amino_acid(&'v'),
            self.get_count_for_amino_acid(&'w'),
            self.get_count_for_amino_acid(&'y')
        ]);
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
            amino_acids_counts: self.amino_acids_counts.clone()
        }
    }
}
