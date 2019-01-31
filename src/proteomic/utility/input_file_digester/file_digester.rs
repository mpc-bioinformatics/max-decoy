pub trait FileDigester {
    fn new(file_path: &str, thread_count: usize, max_number_of_missed_cleavages: i16, min_peptide_length: usize, max_peptide_length: usize) -> Self;
    fn process_file(&mut self, enzym_name: &str) -> f64;
}
