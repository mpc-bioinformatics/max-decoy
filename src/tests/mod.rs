mod models;
mod digest;

#[test]
pub fn run() {
    // peptide tests
    models::peptide::test_equality();
    models::peptide::test_unequlity();
    // protein tests
    models::peptide::test_equality();
    models::peptide::test_unequlity();
    // collection tests
    models::collection::test_len_with_peptides();
    digest::test_digestion();
}