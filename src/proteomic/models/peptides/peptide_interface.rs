use std::collections::HashMap;
use proteomic::models::amino_acids::amino_acid::AMINO_ACIDS_FOR_COUNTING;

/// Trait for accessing fields of a peptide.
/// Peptides are all kinds of amino acids subsequences like Peptide, Decoy, ModifiedPeptides and ModifiedDecoy
pub trait PeptideInterface {
    fn to_string(&self) -> String;
    fn get_aa_sequence(&self) -> &str;
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

    fn as_fasta_entry(header: &str, aa_sequence: &str) -> String {
        return format!("{}\n{}", header, aa_sequence);
    }

    fn count_amino_acids(aa_sequence: &str) -> Box<HashMap<char, i16>> {
        let mut amino_acids_counts: Box<HashMap<char, i16>> = Box::new(HashMap::new());
        for one_letter_code in AMINO_ACIDS_FOR_COUNTING {
            amino_acids_counts.insert(*one_letter_code, aa_sequence.matches(*one_letter_code).count() as i16);
        }
        return amino_acids_counts;
    }
}