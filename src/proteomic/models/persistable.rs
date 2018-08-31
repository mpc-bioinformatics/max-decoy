use super::postgres::Connection;
use super::postgres::rows::Rows;
use super::postgres::Result;
use super::postgres::stmt::Statement;

pub trait Persistable {
    fn get_id(&self) -> i32; 
    fn get_insert_statement() -> &'static str;
    fn execute_insert_query(&self, connection: &Connection) -> Result<Rows>;
    fn execute_insert_statement(&self, prepared_statement: &Statement);
    fn get_update_statement() -> &'static str;
    fn execute_update_query(&self, connection: &Connection) -> Result<Rows>;
    fn execute_update_statement(&self, prepared_statement: &Statement);
    fn get_exists_statement() -> &'static str;
    fn exists_query(&self, connection: &Connection) -> Result<Rows>;
}