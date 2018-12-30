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
}