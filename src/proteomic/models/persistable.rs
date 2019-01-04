extern crate postgres;

use std::string::ToString;
use std::fmt;

pub enum QueryOk {
    Selected,
    Created,
    Updated,
    Deleted,
    AlreadyExists
}

pub enum QueryError {
    NoMatch,
    NoReturn,                           // no return from sql server
    RecordIsNotPersisted,               // if record should updated or delete but is not deleted
    SQLError(String),
    ConnectionError(String),
    IOError(String),
    ConversionError(String),
    InnerError(String),                 // errors which are happend during sql-sresult processing
    AssociatedRecordNotFound(String)
}

impl QueryError {
    fn to_string(&self) -> String {
        return match self {
            QueryError::AssociatedRecordNotFound(err) => format!("QueryError::AssociatedRecordNotFound({})", err),
            QueryError::ConnectionError(err) => format!("QueryError::ConectionError({})", err),
            QueryError::ConversionError(err) => format!("QueryError::ConversionError({})", err),
            QueryError::IOError(err) => format!("QueryError::IOError({})", err),
            QueryError::InnerError(err) => format!("QueryError::InnerError({})", err),
            QueryError::SQLError(err) => format!("QueryError::SqlError(error code: {})", err),
            QueryError::NoMatch => format!("QueryError::NoReturn"),
            QueryError::NoReturn => format!("QueryError::NoMatch"),
            QueryError::RecordIsNotPersisted => format!("QueryError::RecordNotPersisted")
        };
    }
}

impl fmt::Display for QueryError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        return write!(f, "{}", self.to_string());
    }
}

pub enum FromSqlRowError {
    InnerQueryError(QueryError),        // for example if an associated model must be loaded
    AssociatedRecordNotFound(String)
}

impl FromSqlRowError {
    fn to_string(&self) -> String {
        return match self {
            FromSqlRowError::InnerQueryError(err) => format!("FromSqlRowError::InnerQueryErr({})", err),
            FromSqlRowError::AssociatedRecordNotFound(err) => format!("FromSqlRowError::AssociatedRecordNotFound({})", err),
        };
    }
}

impl fmt::Display for FromSqlRowError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        return write!(f, "{}", self.to_string());
    }
}

pub fn handle_postgres_error(error: &postgres::Error) -> QueryError {
    match error.code() {
        Some(sql_state) => return QueryError::SQLError(sql_state.code().to_owned()),
        None => ()
    }
    match error.as_connection() {
        Some(err) => return QueryError::ConnectionError(format!("{}", err)),
        None => ()
    }
    match error.as_conversion() {
        Some(err) => return QueryError::ConversionError(format!("{}", err)),
        None => ()
    }
    match error.as_io() {
        Some(err) => return QueryError::IOError(format!("{}", err)),
        None => ()
    }
    panic!("ERROR [proteomic::models::persistable::handle_postgres_error] some error during query execution happend, which is neither a sql/db-error, nor a connection-error, nor a conversion-erro, nor a io-error")
}

// when implementing
// T is just the type you implement Persistable for
// PK is the type of the primary key
// UI is the type of the unique identifier
pub trait Persistable<T, PK, UI> {
    fn from_sql_row(row: &postgres::rows::Row) -> Result<T, FromSqlRowError>;
    fn get_primary_key(&self) -> PK;
    fn find(conn: &postgres::Connection, primary_key: &PK) -> Result<T, QueryError>;
    fn find_by_unique_identifier(conn: &postgres::Connection, unique_identifier: &UI) -> Result<T, QueryError>;
    fn create(&mut self, conn: &postgres::Connection) -> Result<QueryOk, QueryError>;
    fn update(&mut self, conn: &postgres::Connection) -> Result<QueryOk, QueryError>;
    fn save(&mut self, conn: &postgres::Connection) -> Result<QueryOk, QueryError>;
    fn delete(&mut self, conn: &postgres::Connection) -> Result<QueryOk, QueryError>;
    fn delete_all(conn: &postgres::Connection) -> Result<QueryOk, QueryError>;

    fn select_where(conn: &postgres::Connection, conditions: &str, values: &[&postgres::types::ToSql]) -> Result<Vec<T>, QueryError> {
        let where_statement: String = format!("SELECT * FROM {} WHERE {};", Self::get_table_name(), conditions);
        match conn.query(where_statement.as_str(), values) {
            Ok(ref rows) => {
                let mut records: Vec<T> = Vec::new();
                for row in rows {
                    match Self::from_sql_row(&row) {
                        Ok(record) => records.push(record),
                        Err(err) => return Err(QueryError::InnerError(err.to_string()))
                    }
                    
                }
                return Ok(records);
            },
            Err(err) => Err(handle_postgres_error(&err))
        }
    }

    fn get_table_name() -> &'static str;
    fn get_select_primary_key_by_unique_identifier_query() -> &'static str;
    fn get_insert_query() -> &'static str;
    fn get_update_query() -> &'static str;

    fn exec_select_primary_key_by_unique_identifier_statement(&mut self, prepared_statement: &postgres::stmt::Statement) -> Result<QueryOk, QueryError>;
    fn exec_insert_statement(&mut self, prepared_statement: &postgres::stmt::Statement) -> Result<QueryOk, QueryError>;
    fn exec_update_statement(&mut self, prepared_statement: &postgres::stmt::Statement) -> Result<QueryOk, QueryError>;

    fn is_persisted(&self) -> bool;

    fn get_count(conn: &postgres::Connection) -> Result<i64, QueryError> {
        let query: String = format!("SELECT cast(count(id) AS BIGINT) FROM {};", Self::get_table_name());
        return match conn.query(query.as_str(), &[]) {
            Ok(ref rows) if rows.len() > 0 => Ok(rows.get(0).get::<usize, i64>(0)),
            Ok(_rows) => Err(QueryError::NoReturn),        // this case can not happen, except count is successfull but does not return something. code only exists to make the compiler happy
            Err(err) => Err(handle_postgres_error(&err))
        };
    }
}