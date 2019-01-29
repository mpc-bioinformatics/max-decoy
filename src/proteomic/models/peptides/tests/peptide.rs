use proteomic::models::peptides::peptide::Peptide;

#[test]
pub fn test_equality() {
    let pep1: Peptide = Peptide::new("DHWVHVLVPMGFVIGCYLDR", 0);
    let pep2: Peptide = Peptide::new("DHWVHVLVPMGFVIGCYLDR", 0);
    assert!(pep1 == pep2);
}

#[test]
pub fn test_unequlity() {
    let pep1: Peptide = Peptide::new("DHWVHVLVPMGFVIGCYLDR", 0);
    let pep2: Peptide = Peptide::new("MVNLLQIVR", 0);
    assert!(pep1 != pep2);
}

#[test]
pub fn test_amino_acid_count() {
    let pep: Peptide = Peptide::new("VVGTVK", 0);
    assert_eq!(3, *pep.get_count_for_amino_acid(&'V'));
}