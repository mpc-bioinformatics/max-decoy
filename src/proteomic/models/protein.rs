extern crate onig;
extern crate postgres;

use std::hash::{Hash, Hasher};

use self::postgres::Connection;

use proteomic::models::persistable::{handle_postgres_error, Persistable, QueryError, QueryOk, FromSqlRowError};

pub struct Protein {
    id: i64,                // BIGSERIAL
    accession: String,      // CHAR(10)
    header: String,         // TEXT
    aa_sequence: String     // TEXT
}

impl Protein {
    pub fn new(header: String, aa_sequence: String) -> Protein {
        return Protein {
            id: 0,
            accession: Protein::extract_accession_from_header(&header),
            header: header,
            aa_sequence: aa_sequence
        }
    }

    pub fn is_new(&self) -> bool {
        return self.id < 1;
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

    pub fn as_fasta_entry(&self) -> String {
        return format!("{}\n{}\n", self.header, self.aa_sequence);
    }
}

impl Persistable<Protein, i64, String> for Protein {
    fn from_sql_row(row: &postgres::rows::Row) -> Result<Self, FromSqlRowError> {
        return Ok (
            Self {
                id: row.get(0),
                accession: row.get(1),
                header: row.get(2),
                aa_sequence: row.get(3)
            }
        )
    }

    fn get_primary_key(&self) -> i64 {
        return self.id;
    }

    fn find(conn: &Connection, primary_key: &i64) -> Result<Self, QueryError> {
        match conn.query("SELECT * FROM proteins WHERE id = $1 LIMIT 1", &[primary_key]) {
            Ok(ref rows) if rows.len() > 0 => match Self::from_sql_row(&rows.get(0)) {
                Ok(record) => Ok(record),
                Err(err) => match err {
                    FromSqlRowError::InnerQueryError(from_sql_err) => Err(from_sql_err),
                    FromSqlRowError::AssociatedRecordNotFound(from_sql_err) => Err(QueryError::AssociatedRecordNotFound(from_sql_err.to_string()))
                }
            },
            Ok(_rows) => Err(QueryError::NoMatch),
            Err(err) => Err(handle_postgres_error(&err))
        }
    }

    fn find_by_unique_identifier(conn: &Connection, unique_identifier: &String) -> Result<Self, QueryError> {
        match conn.query("SELECT * FROM proteins WHERE accession = $1 LIMIT 1", &[&unique_identifier]) {
            Ok(ref rows) if rows.len() > 0 => match Self::from_sql_row(&rows.get(0)) {
                Ok(record) => Ok(record),
                Err(err) => match err {
                    FromSqlRowError::InnerQueryError(from_sql_err) => Err(from_sql_err),
                    FromSqlRowError::AssociatedRecordNotFound(from_sql_err) => Err(QueryError::AssociatedRecordNotFound(from_sql_err.to_string()))
                }
            },
            Ok(_rows) => Err(QueryError::NoMatch),
            Err(err) => Err(handle_postgres_error(&err))
        }
    }


    fn create(&mut self, conn: &postgres::Connection) -> Result<QueryOk, QueryError> {
        match conn.query(
            Self::get_insert_query(),
            &[&self.accession, &self.header, &self.aa_sequence]
        ) {
            Ok(ref rows) if rows.len() > 0 => {
                self.id = rows.get(0).get(0);
                return Ok(QueryOk::Created);
            },
            Ok(_rows) => {
                // zero rows means there are a conflict on update, so the protein exists already
                match Self::find_by_unique_identifier(conn, &self.accession) {
                    Ok(protein) => {
                        self.id = protein.get_primary_key();
                        return Ok(QueryOk::AlreadyExists);
                    },
                    Err(err) => Err(err)
                }
            }
            Err(err) => Err(handle_postgres_error(&err))
        }
    }

    fn update(&mut self, conn: &postgres::Connection) -> Result<QueryOk, QueryError> {
        if !self.is_persisted() {
            return Err(QueryError::RecordIsNotPersisted);
        }
        match conn.query(
            Self::get_update_query(),
            &[&self.id, &self.accession, &self.header, &self.aa_sequence]
        ) {
            Ok(ref rows) if rows.len() > 0 => Ok(QueryOk::Updated),
            Ok(_rows) => Err(QueryError::NoReturn),
            Err(err) => Err(handle_postgres_error(&err))
        }
    }

    fn save(&mut self, conn: &postgres::Connection) -> Result<QueryOk, QueryError> {
        if self.is_persisted() {
            return self.update(conn);
        } else {
            return self.create(conn);
        }
    }

    fn delete(&mut self, conn: &postgres::Connection) -> Result<QueryOk, QueryError> {
        if !self.is_persisted() {
            return Err(QueryError::RecordIsNotPersisted);
        }
        match conn.execute("DELETE FROM proteins WHERE id = $1;", &[&self.id]) {
            Ok(_) => {
                self.id = 0;
                return Ok(QueryOk::Deleted);
            },
            Err(err) => Err(handle_postgres_error(&err))
        }
    }

    fn delete_all(conn: &postgres::Connection) -> Result<QueryOk, QueryError> {
        match conn.execute("DELETE FROM proteins WHERE id IS NOT NULL;", &[]) {
            Ok(_) => Ok(QueryOk::Deleted),
            Err(err) => Err(handle_postgres_error(&err))
        }
    }
    

    fn get_table_name() -> &'static str {
        return "proteins";
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


    fn exec_select_primary_key_by_unique_identifier_statement(&mut self, prepared_statement: &postgres::stmt::Statement) -> Result<QueryOk, QueryError> {
        match prepared_statement.query(&[&self.accession]) {
            Ok(ref rows) if rows.len() > 0 => {
                self.id = rows.get(0).get(0);
                return Ok(QueryOk::Selected);
            },
            Ok(_rows) => Err(QueryError::NoMatch),
            Err(err) => Err(handle_postgres_error(&err))
        }
    }

    fn exec_insert_statement(&mut self, prepared_statement: &postgres::stmt::Statement) -> Result<QueryOk, QueryError> {
        match prepared_statement.query(&[&self.accession, &self.header, &self.aa_sequence]) {
            Ok(ref rows) if rows.len() > 0 => {
                self.id = rows.get(0).get(0);
                return Ok(QueryOk::Created);
            },
            Ok(_rows) =>  Err(QueryError::NoReturn),
            Err(err) => Err(handle_postgres_error(&err))
        }
    }

    fn exec_update_statement(&mut self, prepared_statement: &postgres::stmt::Statement) -> Result<QueryOk, QueryError> {
        match prepared_statement.query(&[&self.id, &self.accession, &self.header, &self.aa_sequence]) {
            Ok(ref rows) if rows.len() > 0 => Ok(QueryOk::Updated),
            Ok(_rows) => Err(QueryError::NoReturn),
            Err(err) => Err(handle_postgres_error(&err))
        }
    }


    fn is_persisted(&self) -> bool {
        return self.id > 0;
    }

    fn set_primary_key_from_sql_row(&mut self, row: &postgres::rows::Row) {
        unimplemented!()
    }

    fn exists_query() -> &'static str {
        unimplemented!();
    }

    fn exists_attributes(&self) -> Box<Vec<&postgres::types::ToSql>> {
        unimplemented!();
    }

    fn exists(&mut self, conn: &postgres::Connection) -> Result<QueryOk, QueryError>{
        unimplemented!();
    }

    fn exists_prepared(&mut self, prepared_statement: &postgres::stmt::Statement) -> Result<QueryOk, QueryError>{
        unimplemented!();
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
