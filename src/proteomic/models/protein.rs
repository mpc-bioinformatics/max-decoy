extern crate onig;
extern crate postgres;

use std::hash::{Hash, Hasher};

use self::postgres::Connection;
use self::postgres::rows::Rows;
use self::postgres::Result;
use self::postgres::stmt::Statement;

use proteomic::models::collection::Collectable;
use proteomic::models::persistable::Persistable;

pub struct Protein {
    id: i32,
    accession: String,
    header: String,
    aa_sequence: String,
}

impl Protein {
    pub fn new(header: String, aa_sequence: String) -> Protein {
        return Protein {
            id: -1,
            accession: Protein::extract_accession_from_header(&header),
            header: header,
            aa_sequence: aa_sequence,
        }
    }


    pub fn extract_accession_from_header(header: &String) -> String {
        // return String::from(Regex::new(r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}").unwrap().find(header.as_str()).unwrap().as_str())
        // onigurma has nothing like the original regex::Regex.find
        let accession_regex = onig::Regex::new(r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}").unwrap();
        let pos = accession_regex.find(header.as_str());
        match pos {
            Some((beg, end)) =>
                return String::from(&header[beg..end]),
            None =>
                return String::new()
        }
    }

    pub fn get_aa_sequence(&self) -> &String {
        return &self.aa_sequence;
    }

    pub fn get_accession(&self) -> &String {
        return &self.accession;
    }

    pub fn get_header(&self) -> &String {
        return &self.header;
    }

    pub fn print(&self) {
        println!("{}\n\tAA sequence is {} long.", self.accession, self.aa_sequence.len());
    }

    pub fn create(&mut self, conn: &Connection) {
        for row in self.execute_insert_query(&conn).unwrap().iter() {
            let id: i32 = row.get("id");
            self.id = id;
            break;
        }
    }

    pub fn update(&self, conn: &Connection) {
        self.execute_update_query(&conn);
    }

    pub fn save(&mut self, conn: &Connection) {
        if self.id > 0 {
            self.update(conn);
        } else {
            self.create(conn);
        }
    }

    pub fn exists(&self, conn: &Connection) -> bool {
        for row in self.exists_query(&conn).unwrap().iter() {
            return row.get::<usize, bool>(0);
        }
        return false;
    }
}

impl Persistable for Protein {
    fn get_id(&self) -> i32 {
        return self.id;
    }

    fn get_insert_statement() -> &'static str {
        return "INSERT INTO proteins (accession, header, aa_sequence) VALUES ($1, $2, $3) ON CONFLICT (accession) DO NOTHING";
    }

    fn execute_insert_query(&self, connection: &Connection) -> Result<Rows> {
        return connection.query(
            Protein::get_insert_statement(),
            &[&self.accession, &self.header, &self.aa_sequence]
        );
    }

    fn execute_insert_statement(&self, prepared_statement: &Statement) {
        prepared_statement.execute(&[&self.accession, &self.header, &self.aa_sequence]);
    }

    fn get_update_statement() -> &'static str {
        return "UPDATE proteins SET accession = $2, header = $3, aa_sequence = $4 WHERE id = $1";
    }

    fn execute_update_query(&self, connection: &Connection) -> Result<Rows> {
        return connection.query(
            Protein::get_update_statement(),
            &[&self.id, &self.accession, &self.header, &self.aa_sequence]
        );
    }

    fn execute_update_statement(&self, prepared_statement: &Statement) {
        prepared_statement.execute(&[&self.id, &self.accession, &self.header, &self.aa_sequence]);
    }

    fn get_exists_statement() -> &'static str {
        return "SELECT EXISTS(SELECT 1 FROM proteins WHERE accession = $1)";
    }

    fn exists_query(&self, connection: &Connection) -> Result<Rows> {
        return connection.query(
            Protein::get_exists_statement(),
            &[&self.accession]
        );
    }
}

impl Collectable for Protein {
    fn get_collection_identifier(&self) -> &String {
        return &self.accession;
    }
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