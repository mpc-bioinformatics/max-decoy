extern crate postgres;

// when implementing
// T is just the type you implement Persistable for
// PK is the type of the primary key
// UI is the type of the unique identifier
pub trait Persistable<T, PK, UI> {
    fn from_sql_row(row: &postgres::rows::Row) -> Result<T, String>;
    fn get_primary_key(&self) -> PK;
    fn find(conn: &postgres::Connection, primary_key: &PK) -> Result<T, String>;
    fn find_by_unique_identifier(conn: &postgres::Connection, unique_identifier: &UI) -> Result<T, String>;
    fn create(&mut self, conn: &postgres::Connection) -> Result<(), String>;
    fn update(&mut self, conn: &postgres::Connection) -> Result<(), String>;
    fn save(&mut self, conn: &postgres::Connection) -> Result<(), String>;
    fn delete(&mut self, conn: &postgres::Connection) -> Result<(), String>;
    fn delete_all(conn: &postgres::Connection) -> Result<(), String>;

    fn select_where(conn: &postgres::Connection, conditions: &str, values: &[&postgres::types::ToSql]) -> Result<Vec<T>, String> {
        let where_statement: String = format!("SELECT * FROM {} WHERE {};", Self::get_table_name(), conditions);
        match conn.query(where_statement.as_str(), values) {
            Ok(ref rows) => {
                let mut records: Vec<T> = Vec::new();
                for row in rows {
                    match Self::from_sql_row(&row) {
                        Ok(record) => records.push(record),
                        Err(err) => return Err(err)
                    }
                    
                }
                return Ok(records);
            },
            Err(err) => Err(format!("could not gether records from table '{}'; postgresql error is: {}", Self::get_table_name(), err))
        }
    }

    fn get_table_name() -> &'static str;
    fn get_select_primary_key_by_unique_identifier_query() -> &'static str;
    fn get_insert_query() -> &'static str;
    fn get_update_query() -> &'static str;

    fn exec_select_primary_key_by_unique_identifier_statement(&mut self, prepared_statement: &postgres::stmt::Statement) -> Result<(), String>;
    fn exec_insert_statement(&mut self, prepared_statement: &postgres::stmt::Statement) -> Result<(), String>;
    fn exec_update_statement(&mut self, prepared_statement: &postgres::stmt::Statement) -> Result<(), String>;

    fn is_persisted(&self) -> bool;

    fn get_count(conn: &postgres::Connection) -> Result<i64, String> {
        let query: String = format!("SELECT cast(count(id) AS BIGINT) FROM {};", Self::get_table_name());
        return match conn.query(query.as_str(), &[]) {
            Ok(ref rows) if rows.len() > 0 => Ok(rows.get(0).get::<usize, i64>(0)),
            Ok(_rows) => Err("not possible error occured".to_owned()),                           // this error can not happen, except count is successfull but does not return something. This case exists only to make the compiler happy
            Err(err) => Err(format!("could not count records from table '{}'; original error: {}", Self::get_table_name(), err))
        };
    }
}