mod models;

#[test]
pub fn run() {
    // peptide tests
    models::peptide::test_equality();
    models::peptide::test_unequlity();
    // protein tests
    models::peptide::test_equality();
    models::peptide::test_unequlity();
}