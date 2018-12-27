pub trait Decoy {
    fn to_string(&self) -> String;
    fn get_header(&self) -> String;
    fn get_aa_sequence(&self) -> String;
    fn get_weight(&self) -> i64;
    fn get_length(&self) -> i32;
}