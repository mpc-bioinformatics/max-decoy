use proteomic::models::peptide::Peptide;

#[test]
pub fn test_equality() {
    let pep1: Peptide = Peptide::new(String::from("DHWVHVLVPMGFVIGCYLDR"), String::from("Trypsin"), 0);
    let pep2: Peptide = Peptide::new(String::from("DHWVHVLVPMGFVIGCYLDR"), String::from("Trypsin"), 0);
    assert!(pep1 == pep2);
}

#[test]
pub fn test_unequlity() {
    let pep1: Peptide = Peptide::new(String::from("DHWVHVLVPMGFVIGCYLDR"), String::from("Trypsin"), 0);
    let pep2: Peptide = Peptide::new(String::from("MVNLLQIVR"), String::from("Trypsin"), 0);
    assert!(pep1 != pep2);
}