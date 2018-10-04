extern crate postgres;

pub trait Persistable<T, PK, UI> {
    fn get_primary_key(&self) -> &PK;
    fn find(conn: &postgres::Connection, primary_key: &PK) -> Result<T, &'static str>;
    fn find_by_unique_identifier(conn: &postgres::Connection, unique_identifier: &UI) -> Result<T, &'static str>;
    fn create(&self, conn: &postgres::Connection) -> bool;
    fn update(&self, conn: &postgres::Connection) -> bool;
    fn save(&self, conn: &postgres::Connection) -> bool;

    fn get_insert_query() -> &'static str;
    fn get_update_query() -> &'static str;

    fn exec_insert_statement(&self, prepared_statement: &postgres::stmt::Statement) -> bool;
    fn exec_update_statement(&self, prepared_statement: &postgres::stmt::Statement) -> bool;

    fn is_persisted(&self) -> bool;
}