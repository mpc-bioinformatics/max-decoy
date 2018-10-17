extern crate dotenv;
extern crate postgres;

use std::env;

pub struct DatabaseConnection {}

impl DatabaseConnection {
    fn get_database_url() -> String {
        dotenv::dotenv().ok();
        return env::var("PGSQL_URL").expect("Postgresql-Database-URL 'PGSQL_URL' must be set ");
    }

    pub fn get_database_connection() -> postgres::Connection {
        return postgres::Connection::connect(DatabaseConnection::get_database_url().as_str(), postgres::TlsMode::None).unwrap();
    }
}