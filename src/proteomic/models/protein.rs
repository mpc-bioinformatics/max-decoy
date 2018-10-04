extern crate onig;
extern crate postgres;

use std::hash::{Hash, Hasher};

use self::postgres::Connection;
use self::postgres::rows::Rows;
use self::postgres::stmt::Statement;

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
    fn get_primary_key(&self) -> &i32 {
        return &self.id;
    }

    fn find(conn: &Connection, primary_key: &i32) -> Result<Self, &'static str> {
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
                    Err("peptide not found")
                }
            },
            Err(_err) => Err("some error occured: Peptide::find")
        }
    }

    fn find_by_unique_identifier(conn: &Connection, unique_identifier: &String) -> Result<Self, &'static str> {
        match conn.query("SELECT * FROM proteins WHERE aa_sequence = $1 LIMIT 1", &[&unique_identifier]) {
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
                    Err("peptide not found")
                }
            },
            Err(_err) => Err("some error occured: Peptide::find_by_unique_identifier")
        }
    }


    fn create(&self, conn: &postgres::Connection) -> bool {
        match conn.query(
            Self::get_insert_query(),
            &[&self.accession, &self.header, &self.aa_sequence]
        ) {
            Ok(rows) => {
                if rows.len() > 0 {
                    self.id =  rows.get(0).get(0);
                    return true;
                } else {
                    // zero rows means there are a conflict on update, so the peptides exists already
                    match Self::find_by_unique_identifier(conn, &self.aa_sequence) {
                        Ok(peptide) => {
                            self.id = *peptide.get_primary_key();
                            self.is_persisted = true;
                            return true;
                        },
                        Err(_err) => {
                            println!("ERROR: something is really odd with peptide {}:\n\tcan't create it nor find it\n", self.aa_sequence);
                            return false;
                        }
                    }
                }
            }
            Err(_err) => return false
        }
    }

    fn update(&self, conn: &postgres::Connection) -> bool {
        match conn.query(
            Self::get_update_query(),
            &[&self.id, &self.accession, &self.header, &self.aa_sequence]
        ) {
            Ok(rows) => return rows.len() > 0,
            Err(_err) => return false
        }
    }

    fn save(&self, conn: &postgres::Connection) -> bool {
        if !self.is_persisted() {
            return self.update(conn);
        } else {
            return self.create(conn);
        }
    }


    fn get_insert_query() -> &'static str {
        return "INSERT INTO proteins (accession, header, aa_sequence) VALUES ($1, $2, $3) ON CONFLICT DO NOTHING RETURNING id";
    }

    fn get_update_query() -> &'static str {
        return "UPDATE proteins SET accession = $2, header = $3, aa_sequence = $4 WHERE id = $1";
    }


    fn exec_insert_statement(&self, prepared_statement: &postgres::stmt::Statement) -> bool {
        match prepared_statement.query(&[&self.accession, &self.header, &self.aa_sequence]) {
            Ok(rows) => {
                if rows.len() > 0 {
                    self.id = rows.get(0).get(0);
                    self.is_persisted = true;
                    return true;
                } else {
                    return false
                }
            },
            Err(_err) => return false
        }
    }

    fn exec_update_statement(&self, prepared_statement: &postgres::stmt::Statement) -> bool {
        match prepared_statement.query(&[&self.id, &self.accession, &self.header, &self.aa_sequence]) {
            Ok(rows) => return rows.len() > 0,
            Err(_err) => return false
        }
    }


    fn is_persisted(&self) -> bool {
        return self.is_persisted;
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