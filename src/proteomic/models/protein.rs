extern crate regex;

pub struct Protein{
    accession: String,
    header: String,
    aa_sequence: String
}

impl Protein{
    pub fn new(header: String, aa_sequence: String) -> Protein {
        return Protein{
            accession: Protein::extract_accession_from_header(&header),
            header: header,
            aa_sequence: aa_sequence
        }
    }

    pub fn extract_accession_from_header(header: &String) -> String {
        return String::from(regex::Regex::new(r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}").unwrap().find(header.as_str()).unwrap().as_str())
    }

    pub fn get_aa_sequence(&self) -> &String {
        &self.aa_sequence
    }

    pub fn get_accession(&self) -> &String {
        &self.accession
    }

    pub fn get_header(&self) -> &String {
        &self.accession
    }

    pub fn print(&self) {
        println!("{}\n\tAA sequence is {} long.", self.accession, self.aa_sequence.len());
    }
}