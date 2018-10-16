extern crate onig;
extern crate postgres;

use std::hash::{Hash, Hasher};
use std::error::Error;

use self::postgres::Connection;

use proteomic::models::collection::Collectable;
use proteomic::models::persistable::Persistable;

pub struct Protein {
    id: i32,
    accession: String,
    header: String,
    aa_sequence: String,
    is_persisted: bool
}

impl Protein {
    pub fn new(header: String, aa_sequence: String) -> Protein {
        return Protein {
            id: -1,
            accession: Protein::extract_accession_from_header(&header),
            header: header,
            aa_sequence: aa_sequence,
            is_persisted: false
        }
    }

    pub fn is_new(&self) -> bool {
        return self.id < 0;
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

    pub fn to_string(&self) -> String {
        return format!("{}: {}\n\tlen => {}", self.id, self.accession, self.aa_sequence.len());
    }
}

impl Persistable<Protein, i32, String> for Protein {
    fn get_primary_key(&self) -> i32 {
        return self.id;
    }

    fn find(conn: &Connection, primary_key: &i32) -> Result<Self, String> {
        match conn.query("SELECT * FROM proteins WHERE id = $1 LIMIT 1", &[primary_key]) {
            Ok(rows) =>{
                if rows.len() > 0 {
                    Ok(
                        Protein{
                            id: rows.get(0).get(0),
                            accession: rows.get(0).get(1),
                            header: rows.get(0).get(2),
                            aa_sequence: rows.get(0).get(3),
                            is_persisted: true
                        }
                    )
                } else {
                    Err("NOHIT".to_owned())
                }
            },
            Err(err) => Err(err.code().unwrap().code().to_owned())
        }
    }

    fn find_by_unique_identifier(conn: &Connection, unique_identifier: &String) -> Result<Self, String> {
        match conn.query(
            "SELECT * FROM proteins WHERE accession = $1 LIMIT 1",
            &[&unique_identifier]
        ) {
            Ok(ref rows)  if rows.len() > 0 =>Ok(
                Protein{
                    id: rows.get(0).get(0),
                    accession: rows.get(0).get(1),
                    header: rows.get(0).get(2),
                    aa_sequence: rows.get(0).get(3),
                    is_persisted: true
                }
            ),
            Ok(_rows) => Err("NOHIT".to_owned()),
            Err(err) => Err(err.code().unwrap().code().to_owned())
        }
    }


    fn create(&mut self, conn: &postgres::Connection) -> Result<(), String> {
        match conn.query(
            Self::get_insert_query(),
            &[&self.accession, &self.header, &self.aa_sequence]
        ) {
            Ok(ref rows) if rows.len() > 0 => {
                self.id = rows.get(0).get(0);
                return Ok(());
            },
            Ok(_rows) => {
                // zero rows means there are a conflict on update, so the protein exists already
                match Self::find_by_unique_identifier(conn, &self.accession) {
                    Ok(protein) => {
                        self.id = protein.get_primary_key();
                        self.is_persisted = true;
                        return Ok(());
                    },
                    Err(err) => Err(format!("cannot insert nor find protein '{}'\n\torginal error {}", self.accession, err))
                }
            }
            Err(err) => Err(err.code().unwrap().code().to_owned())
        }
    }

    fn update(&mut self, conn: &postgres::Connection) -> Result<(), String> {
        match conn.query(
            Self::get_update_query(),
            &[&self.id, &self.accession, &self.header, &self.aa_sequence]
        ) {
            Ok(ref rows) if rows.len() > 0 => Ok(()),
            Ok(_rows) => Err("NORET".to_owned()),
            Err(err) => Err(err.code().unwrap().code().to_owned())
        }
    }

    fn save(&mut self, conn: &postgres::Connection) -> Result<(), String> {
        if self.is_persisted() {
            return self.update(conn);
        } else {
            return self.create(conn);
        }
    }


    fn get_select_primary_key_by_unique_identifier_query() -> &'static str {
        return "SELECT id FROM proteins WHERE accession = $1 LIMIT 1";
    }

    fn get_insert_query() -> &'static str {
        return "INSERT INTO proteins (accession, header, aa_sequence) VALUES ($1, $2, $3) ON CONFLICT DO NOTHING RETURNING id";
    }

    fn get_update_query() -> &'static str {
        return "UPDATE proteins SET accession = $2, header = $3, aa_sequence = $4 WHERE id = $1";
    }


    fn exec_select_primary_key_by_unique_identifier_statement(&mut self, prepared_statement: &postgres::stmt::Statement) -> Result<(), String> {
        match prepared_statement.query(&[&self.accession]) {
            Ok(ref rows) if rows.len() > 0 => {
                self.id = rows.get(0).get(0);
                self.is_persisted = true;
                return Ok(());
            },
            Ok(_rows) => Err("NOHIT".to_owned()),
            Err(err) => Err(err.code().unwrap().code().to_owned())
        }
    }

    fn exec_insert_statement(&mut self, prepared_statement: &postgres::stmt::Statement) -> Result<(), String> {
        match prepared_statement.query(&[&self.accession, &self.header, &self.aa_sequence]) {
            Ok(ref rows) if rows.len() > 0 => {
                self.id = rows.get(0).get(0);
                self.is_persisted = true;
                return Ok(());
            },
            Ok(_rows) =>  Err("NORET".to_owned()),
            Err(err) => Err(err.code().unwrap().code().to_owned())
        }
    }

    fn exec_update_statement(&mut self, prepared_statement: &postgres::stmt::Statement) -> Result<(), String> {
        match prepared_statement.query(&[&self.id, &self.accession, &self.header, &self.aa_sequence]) {
            Ok(ref rows) if rows.len() > 0 => Ok(()),
            Ok(_rows) => Err("NORET".to_owned()),
            Err(err) => Err(err.code().unwrap().code().to_owned())
        }
    }


    fn is_persisted(&self) -> bool {
        return self.is_persisted;
    }

    fn get_count(conn: &postgres::Connection) -> i64 {
        return match conn.query("SELECT cast(count(id) AS BIGINT) FROM proteins", &[]) {
            Ok(ref rows) if rows.len() > 0 => rows.get(0).get::<usize, i64>(0),
            _ => -1
        };
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