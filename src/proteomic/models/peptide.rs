use std::collections::HashSet;
use std::hash::{Hash, Hasher};

pub struct Peptide {
    id: i64,
    protein_ids: HashSet<i64>,
    aa_sequence: String,
    position_in_protein: usize,
    digest_enzym: String,
    length: usize,
}

impl Peptide {
    pub fn new(aa_sequence: String, position_in_protein: usize, digest_enzym: String) -> Peptide {
        return Peptide{
            id: -1,
            protein_ids: HashSet::new(),
            length: aa_sequence.len(),
            aa_sequence: aa_sequence,
            position_in_protein: position_in_protein,
            digest_enzym: digest_enzym,
        }
    }

    pub fn print(&self) {
        println!("{} => digested with {} : found in #{}", self.aa_sequence, self.digest_enzym, self.protein_ids.len());
    }
}

// PartialEq-implementation to use this type in a HashSet
impl PartialEq for Peptide {
    fn eq(&self, other: &Peptide) -> bool {
       return self.aa_sequence.eq(&other.aa_sequence);
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