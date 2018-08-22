use std::hash::{Hash, Hasher};

pub struct Peptide {
    aa_sequence: String,
    digest_enzym: String,
    length: usize,
    position: usize,
    number_of_missed_cleavages: usize
}

impl Peptide {
    pub fn new(aa_sequence: String, digest_enzym: String, position: usize, number_of_missed_cleavages: usize) -> Peptide {
        return Peptide{
            length: aa_sequence.len(),
            aa_sequence: aa_sequence,
            digest_enzym: digest_enzym,
            position: position,
            number_of_missed_cleavages: number_of_missed_cleavages
        }
    }

    pub fn print(&self) {
        println!("{}\n\tdigested with => {}\n\tpos => {}\n\tmissed_cleavages => {}\n", self.aa_sequence, self.digest_enzym, self.position, self.number_of_missed_cleavages);
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