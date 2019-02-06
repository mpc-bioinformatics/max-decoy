use std::hash::{Hash, Hasher};

pub struct FastaEntry {
    header: String,
    aa_sequence: String
}

impl FastaEntry {
    pub fn new(header: &str, aa_sequence: &str) -> Self {
        return Self {
            header: header.to_owned(),
            aa_sequence: aa_sequence.to_owned()
        }
    }

    pub fn to_string(&self) -> String {
        return format!("{}\n{}\n", self.header.as_str(), self.aa_sequence.as_str())
    }

    pub fn get_header(&self) -> &str {
        return self.header.as_str();
    }

    pub fn get_aa_sequence(&self) -> &str {
        return self.aa_sequence.as_str();
    }
}

// PartialEq-implementation to use this type in a HashSet
impl PartialEq for FastaEntry {
    fn eq(&self, other: &FastaEntry) -> bool {
        return self.aa_sequence.eq(other.get_aa_sequence());
    }
}

// Eq-implementation to use this type in a HashSet
impl Eq for FastaEntry {}

// Hash-implementation to use this type in a HashSet
impl Hash for FastaEntry {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.aa_sequence.hash(state);
    }
}

impl Clone for FastaEntry {
    fn clone(&self) -> FastaEntry {
        return FastaEntry {
            header: self.header.clone(),
            aa_sequence: self.aa_sequence.clone()
        }
    }
}