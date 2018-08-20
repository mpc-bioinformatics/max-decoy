pub struct Protein{
    accession: String,
    aa_sequence: String
}

impl Protein{
    pub fn new(accession: String, aa_sequence: String) -> Protein {
        Protein{
            accession: accession,
            aa_sequence: aa_sequence
        }
    }

    pub fn get_aa_sequence(&self) -> &String {
        &self.aa_sequence
    }

    pub fn get_accession(&self) -> &String {
        &self.accession
    }

    pub fn print(&self) {
        println!("{}\nIt's aa sequence is {} long.", self.accession, self.aa_sequence.len());
    }
}