
use std::string::ToString;

/// Trait for accessing fields of a peptide.
/// Peptides are all kinds of amino acids subsequences like Peptide, Decoy, ModifiedPeptides and ModifiedDecoy
pub trait PeptideInterface {
    fn to_string(&self) -> String;
    fn get_header(&self) -> &str;
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