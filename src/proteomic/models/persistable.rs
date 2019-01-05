extern crate postgres;

use std::fmt;

pub enum QueryOk {
    Selected,
    Created,
    Updated,
    Deleted,
    Exists,
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
            QueryError::NoMatch => format!("QueryError::NoMatch"),
            QueryError::NoReturn => format!("QueryError::NoReturn"),
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
    ////// utilities
    fn from_sql_row(row: &postgres::rows::Row) -> Result<T, FromSqlRowError>;
    fn set_primary_key_from_sql_row(&mut self, row: &postgres::rows::Row);
    fn invalidate_primary_key(&mut self);
    fn get_primary_key(&self) -> PK;
    fn get_table_name() -> &'static str;
    fn is_persisted(&self) -> bool;
    /////

    /////// queries
    /// find-query
    fn find_query() -> &'static str;
    fn find(conn: &postgres::Connection, primary_key: &[&postgres::types::ToSql]) -> Result<T, QueryError> {
        match conn.query(Self::find_query(), primary_key) {
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
    ///
    
    /// create-query
    fn create_query() -> &'static str;
    fn create_attributes(&self) -> Box<Vec<&postgres::types::ToSql>>;
    fn create(&mut self, conn: &postgres::Connection) -> Result<QueryOk, QueryError> {
        match self.exists(&conn) {
            Ok(_) => return Ok(QueryOk::AlreadyExists),
            Err(err) => match err {
                QueryError::NoMatch => match conn.query( Self::create_query(), self.create_attributes().as_slice()) {
                    Ok(ref rows) if rows.len() > 0 => {
                        self.set_primary_key_from_sql_row(&rows.get(0));
                        return Ok(QueryOk::Created);
                    },
                    Ok(_rows) => Err(QueryError::NoReturn),
                    Err(err) => Err(handle_postgres_error(&err))
                },
                _ => Err(err)
            }
        }
    }
    fn prepared_create(&mut self, create_statement: &postgres::stmt::Statement, exists_statement: &postgres::stmt::Statement) -> Result<QueryOk, QueryError> {
        match self.prepared_exists(exists_statement) {
            Ok(_) => return Ok(QueryOk::AlreadyExists),
            Err(err) => match err {
                QueryError::NoMatch => match create_statement.query(self.create_attributes().as_slice()) {
                    Ok(ref rows) if rows.len() > 0 => {
                        self.set_primary_key_from_sql_row(&rows.get(0));
                        return Ok(QueryOk::Created);
                    },
                    Ok(_rows) => Err(QueryError::NoReturn),
                    Err(err) => Err(handle_postgres_error(&err))
                },
                _ => Err(err)
            }
        }
    }
    ///
    
    /// update
    fn update_query() -> &'static str;
    fn update_attributes(&self) -> Box<Vec<&postgres::types::ToSql>>;
    fn update(&mut self, conn: &postgres::Connection) -> Result<QueryOk, QueryError>{
        if !self.is_persisted() {
            return Err(QueryError::RecordIsNotPersisted);
        }
        match conn.query(Self::update_query(), self.update_attributes().as_slice()) {
            Ok(ref rows) if rows.len() > 0 => Ok(QueryOk::Updated),
            Ok(_rows) => Err(QueryError::NoReturn),
            Err(err) => Err(handle_postgres_error(&err))
        }
    }
    fn prepared_update(&mut self, prepared_statement: &postgres::stmt::Statement) -> Result<QueryOk, QueryError> {
        match prepared_statement.query(self.update_attributes().as_slice()) {
            Ok(ref rows) if rows.len() > 0 => Ok(QueryOk::Updated),
            Ok(_rows) => Err(QueryError::NoReturn),
            Err(err) => Err(handle_postgres_error(&err))
        }
    }
    ///
    
    
    /// delete
    fn delete_query() -> &'static str;
    fn delete_attributes(&self) -> Box<Vec<&postgres::types::ToSql>>;
    fn delete(&mut self, conn: &postgres::Connection) -> Result<QueryOk, QueryError> {
        if !self.is_persisted() {
            return Err(QueryError::RecordIsNotPersisted);
        }
        match conn.execute(Self::delete_query(), self.delete_attributes().as_slice()) {
            Ok(_) => {
                self.invalidate_primary_key();
                return Ok(QueryOk::Deleted);
            },
            Err(err) => Err(handle_postgres_error(&err))
        }
    }
    ///

    /// delete_all
    fn delete_all_query() -> &'static str;
    fn delete_all(conn: &postgres::Connection) -> Result<QueryOk, QueryError> {
        match conn.execute(Self::delete_all_query(), &[]) {
            Ok(_) => Ok(QueryOk::Deleted),
            Err(err) => Err(handle_postgres_error(&err))
        }
    }
    ///
    
    /// delete where
    fn delete_where(conn: &postgres::Connection, conditions: &str, values: &[&postgres::types::ToSql]) -> Result<QueryOk, QueryError> {
        let delete_query: String = format!("DELETE FROM {} WHERE {};", Self::get_table_name(), conditions);
        match conn.execute(delete_query.as_str(), values) {
            Ok(_) => Ok(QueryOk::Deleted),
            Err(err) => Err(handle_postgres_error(&err))
        }
    }
    /// 
 
    /// exists-query
    fn exists_query() -> &'static str;
    fn exists_attributes(&self) -> Box<Vec<&postgres::types::ToSql>>;
    fn exists(&mut self, conn: &postgres::Connection) -> Result<QueryOk, QueryError> {
        match conn.query(Self::exists_query(), self.exists_attributes().as_slice()) {
            Ok(ref rows) if rows.len() > 0 => {
                self.set_primary_key_from_sql_row(&rows.get(0));
                return Ok(QueryOk::Exists);
            },
            Ok(_rows) => return Err(QueryError::NoMatch),
            Err(err) => Err(handle_postgres_error(&err))
        }
    }
    fn prepared_exists(&mut self, prepared_statement: &postgres::stmt::Statement) -> Result<QueryOk, QueryError> {
        match prepared_statement.query(self.exists_attributes().as_slice()) {
            Ok(ref rows) if rows.len() > 0 => {
                self.set_primary_key_from_sql_row(&rows.get(0));
                return Ok(QueryOk::Exists);
            },
            Ok(_rows) => return Err(QueryError::NoMatch),
            Err(err) => Err(handle_postgres_error(&err))
        }
    }
    ///
    
    /// select 
    fn select_where(conn: &postgres::Connection, conditions: &str, values: &[&postgres::types::ToSql]) -> Result<Vec<T>, QueryError> {
        let select_query: String = format!("SELECT * FROM {} WHERE {};", Self::get_table_name(), conditions);
        match conn.query(select_query.as_str(), values) {
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
    ///
    
    /// save
    fn save(&mut self, conn: &postgres::Connection) -> Result<QueryOk, QueryError> {
        if self.is_persisted() {
            return self.update(conn);
        } else {
            return self.create(conn);
        }
    }
    /// 
    
    /// count 
    fn get_count(conn: &postgres::Connection) -> Result<i64, QueryError> {
        let query: String = format!("SELECT cast(count(id) AS BIGINT) FROM {};", Self::get_table_name());
        return match conn.query(query.as_str(), &[]) {
            Ok(ref rows) if rows.len() > 0 => Ok(rows.get(0).get::<usize, i64>(0)),
            Ok(_rows) => Err(QueryError::NoReturn),
            Err(err) => Err(handle_postgres_error(&err))
        };
    }
    ///
    //////
    
    ////// hooks
    /// 
    fn before_delete_hook(&self) -> Result<(), QueryError>;
    //
}