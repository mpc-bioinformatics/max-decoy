extern crate postgres;

pub trait Persistable<  T, PK, UI> {
    fn get_primary_key(&self) -> PK;
    fn find(conn: &postgres::Connection, primary_key: &PK) -> Result<T, String>;
    fn find_by_unique_identifier(conn: &postgres::Connection, unique_identifier: &UI) -> Result<T, String>;
    fn create(&mut self, conn: &postgres::Connection) -> Result<(), String>;
    fn update(&mut self, conn: &postgres::Connection) -> Result<(), String>;
    fn save(&mut self, conn: &postgres::Connection) -> Result<(), String>;

    fn get_select_primary_key_by_unique_identifier_query() -> &'static str;
    fn get_insert_query() -> &'static str;
    fn get_update_query() -> &'static str;

    fn exec_select_primary_key_by_unique_identifier_statement(&mut self, prepared_statement: &postgres::stmt::Statement) -> Result<(), String>;
    fn exec_insert_statement(&mut self, prepared_statement: &postgres::stmt::Statement) -> Result<(), String>;
    fn exec_update_statement(&mut self, prepared_statement: &postgres::stmt::Statement) -> Result<(), String>;

    fn is_persisted(&self) -> bool;

    fn get_count(conn: &postgres::Connection) -> i64;
}