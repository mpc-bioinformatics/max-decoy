extern crate postgres;

use self::postgres::Connection;

use proteomic::models::persistable::{handle_postgres_error, Persistable, QueryError, QueryOk, FromSqlRowError};
use proteomic::models::peptide::Peptide;
use proteomic::models::protein::Protein;

pub struct PeptideProteinAssociation {
    peptide_id: i64,
    protein_id: i64
}

impl PeptideProteinAssociation {
    pub fn new(peptide: &Peptide, protein: &Protein) -> PeptideProteinAssociation {
        return PeptideProteinAssociation {
            peptide_id: peptide.get_primary_key(),
            protein_id: protein.get_primary_key()
        }
    }

    pub fn get_peptide_id(&self) -> i64 {
        return self.peptide_id;
    }

    pub fn get_protein_id(&self) -> i64 {
        return self.protein_id;
    }
}

impl Persistable<PeptideProteinAssociation, (i64, i64), (i64, i64)> for PeptideProteinAssociation {
    fn from_sql_row(row: &postgres::rows::Row) -> Result<Self, FromSqlRowError> {
        return Ok(
            Self {
                peptide_id: row.get(0),
                protein_id: row.get(1)
            }
        )
    }

    fn get_primary_key(&self) -> (i64, i64) {
        return (self.peptide_id, self.protein_id);
    }

    fn find(conn: &Connection, primary_key: &(i64, i64)) -> Result<Self, QueryError> {
        return Self::find_by_unique_identifier(conn, primary_key);
    }

    fn find_by_unique_identifier(conn: &Connection, unique_identifier: &(i64, i64)) -> Result<Self, QueryError> {
        match conn.query("SELECT * FROM peptides WHERE peptide_id = $1 and protein_id = $2 LIMIT 1", &[&unique_identifier.0, &unique_identifier.1]) {
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
            &[&self.peptide_id, &self.protein_id]
        ) {
            Ok(ref rows) if rows.len() > 0 => {
                self.peptide_id = rows.get(0).get(0);
                self.protein_id = rows.get(0).get(1);
                return Ok(QueryOk::Created);
            }
            Ok(_rows) => {
                // zero rows means there are a conflict on update, so the peptides exists already
                match Self::find_by_unique_identifier(conn, &(self.peptide_id, self.protein_id)) {
                    Ok(peptide_protein_association) => {
                        self.peptide_id = peptide_protein_association.get_primary_key().0;
                        self.protein_id = peptide_protein_association.get_primary_key().1;
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
            &[&self.peptide_id, &self.protein_id]
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
        match conn.execute("DELETE FROM peptides_proteins WHERE peptide_id = $1 AND protein_id = $2;", &[&self.peptide_id, &self.protein_id]) {
            Ok(_) => {
                self.peptide_id = 0;
                self.protein_id = 0;
                return Ok(QueryOk::Deleted);
            },
            Err(err) => Err(handle_postgres_error(&err))
        }
    }

    fn delete_all(conn: &postgres::Connection) -> Result<QueryOk, QueryError> {
        match conn.execute("DELETE FROM peptides_proteins WHERE peptide_id IS NOT NULL AND protein_id IS NOT NULL;", &[]) {
            Ok(_) => Ok(QueryOk::Deleted),
            Err(err) => Err(handle_postgres_error(&err))
        }
    }


    fn get_table_name() -> &'static str {
        return "peptides_proteins";
    }

    fn get_select_primary_key_by_unique_identifier_query() -> &'static str {
        return "SELECT peptide_id, protein_id FROM peptides_proteins WHERE peptide_id = $1 and protein_id = $2 LIMIT 1";
    }

    fn get_insert_query() -> &'static str {
        return "INSERT INTO peptides_proteins (peptide_id, protein_id) VALUES ($1, $2) ON CONFLICT (peptide_id, protein_id) DO NOTHING RETURNING *";
    }

    fn get_update_query() -> &'static str {
        return "UPDATE peptides_proteins SET peptide_id = $1 AND protein_id = $2 WHERE peptide_id = $1 AND protein_id = $2";
    }


    fn exec_select_primary_key_by_unique_identifier_statement(&mut self, prepared_statement: &postgres::stmt::Statement) -> Result<QueryOk, QueryError> {
        match prepared_statement.query(&[&self.peptide_id, &self.protein_id]) {
            Ok(ref rows) if rows.len() > 0 => {
                self.peptide_id = rows.get(0).get(0);
                self.protein_id = rows.get(0).get(1);
                return Ok(QueryOk::Selected);
            },
            Ok(_rows) => Err(QueryError::NoReturn),
            Err(err) => Err(handle_postgres_error(&err))
        }
    }

    fn exec_insert_statement(&mut self, prepared_statement: &postgres::stmt::Statement) -> Result<QueryOk, QueryError> {
        match prepared_statement.query(&[&self.peptide_id, &self.protein_id]) {
            Ok(ref rows) if rows.len() > 0 => {
                self.peptide_id = rows.get(0).get(0);
                self.protein_id = rows.get(0).get(1);
                return Ok(QueryOk::Created);
            },
            Ok(_rows) => Err(QueryError::NoReturn),
            Err(err) => Err(handle_postgres_error(&err))
        }
    }

    fn exec_update_statement(&mut self, prepared_statement: &postgres::stmt::Statement) -> Result<QueryOk, QueryError> {
        match prepared_statement.query(&[&self.peptide_id, &self.protein_id]) {
            Ok(ref rows) if rows.len() > 0 => Ok(QueryOk::Updated),
            Ok(_rows) => Err(QueryError::NoReturn),
            Err(err) => Err(handle_postgres_error(&err))
        }
    }


    fn is_persisted(&self) -> bool {
        return (self.peptide_id > 0) & (self.protein_id > 0);
    }
}
