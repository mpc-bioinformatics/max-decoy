use std::hash::{Hash, Hasher};

use proteomic::models::persistable::{Persistable, QueryError, FromSqlRowError};

pub struct Protein {
    id: i64,                        // BIGSERIAL
    accession: String,              // CHAR(10)
    header: String,                 // TEXT
    aa_sequence: String,            // TEXT
    is_completely_digested: bool    // BOOLEAN
}

impl Protein {
    pub fn new(header: &str, aa_sequence: &str) -> Protein {
        return Protein {
            id: 0,
            accession: Protein::extract_accession_from_header(header),
            header: header.to_owned(),
            aa_sequence: aa_sequence.to_owned(),
            is_completely_digested: false
        }
    }

    pub fn extract_accession_from_header(header: &str) -> String {
        // return String::from(Regex::new(r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}").unwrap().find(header.as_str()).unwrap().as_str())
        // onigurma has nothing like the original regex::Regex.find
        let accession_regex = onig::Regex::new(r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}").unwrap();
        let pos = accession_regex.find(header);
        match pos {
            Some((beg, end)) =>
                return String::from(&header[beg..end]),
            None =>
                return String::new()
        }
    }

    pub fn get_aa_sequence(&self) -> &str {
        return self.aa_sequence.as_str();
    }

    pub fn get_accession(&self) -> &str {
        return self.accession.as_str();
    }

    pub fn get_header(&self) -> &str {
        return self.header.as_str();
    }

    pub fn to_string(&self) -> String {
        return format!("{}: {}\n\tlen => {}", self.id, self.accession, self.aa_sequence.len());
    }

    pub fn as_fasta_entry(&self) -> String {
        return format!("{}\n{}\n", self.header, self.aa_sequence);
    }

    pub fn is_completely_digested(&self) -> bool {
        return self.is_completely_digested;
    }

    pub fn set_is_completely_digested(&mut self, is_completely_digested: bool) {
        self.is_completely_digested = is_completely_digested;
    }
}

impl Persistable<Protein, i64, String> for Protein {
    fn from_sql_row(row: &postgres::rows::Row) -> Result<Self, FromSqlRowError> {
        return Ok (
            Self {
                id: row.get(0),
                accession: row.get(1),
                header: row.get(2),
                aa_sequence: row.get(3),
                is_completely_digested: row.get(4)
            }
        )
    }

    fn set_primary_key_from_sql_row(&mut self, row: &postgres::rows::Row) {
        self.id = row.get(0);
    }

    fn invalidate_primary_key(&mut self) {
        self.id = 0;
    }

    fn get_primary_key(&self) -> i64 {
        return self.id;
    }

    fn get_table_name() -> &'static str {
        return "proteins";
    }

    fn is_persisted(&self) -> bool {
        return self.id > 0;
    }

    fn find_query() -> &'static str {
        return "SELECT * FROM proteins WHERE id = $1 LIMIT 1;";
    }

    fn create_query() -> &'static str {
        return "INSERT INTO proteins (accession, header, aa_sequence, is_completely_digested) VALUES ($1, $2, $3, $4) ON CONFLICT DO NOTHING RETURNING id;";
    }

    fn create_attributes(&self) -> Box<Vec<&postgres::types::ToSql>>{
        return Box::new(vec![&self.accession, &self.header, &self.aa_sequence, &self.is_completely_digested]);
    }

    fn update_query() -> &'static str{
        return "UPDATE proteins SET accession = $2, header = $3, aa_sequence = $4, is_completely_digested = $5 WHERE id = $1;";
    }

    fn update_attributes(&self) -> Box<Vec<&postgres::types::ToSql>>{
        return Box::new(vec![&self.id, &self.accession, &self.header, &self.aa_sequence, &self.is_completely_digested]);
    }

    fn delete_query() -> &'static str {
        return "DELETE FROM proteins WHERE id = $1;";
    }

    fn delete_attributes(&self) -> Box<Vec<&postgres::types::ToSql>> {
        return Box::new(vec![&self.id]);
    }

    fn delete_all_query() -> &'static str {
        return "DELETE FROM proteins WHERE id IS NOT NULL;";
    }

    fn exists_query() -> &'static str {
        return "SELECT id FROM proteins WHERE accession = $1 LIMIT 1;";
    }

    fn exists_attributes(&self) -> Box<Vec<&postgres::types::ToSql>> {
        return Box::new(vec![&self.accession]);
    }

    fn before_delete_hook(&self) -> Result<(), QueryError> {return Ok(());}
}

// PartialEq-implementation to use this type in a HashSet
impl PartialEq for Protein {
    fn eq(&self, other: &Protein) -> bool {
       return self.accession.eq(&other.accession);
    }
}

// Eq-implementation to use this type in a HashSet
impl Eq for Protein {}

// Hash-implementation to use this type in a HashSet
impl Hash for Protein {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.accession.hash(state);
    }
}
