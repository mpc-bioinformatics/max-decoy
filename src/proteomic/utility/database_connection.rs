use std::env;

use proteomic::models::persistable::handle_postgres_error;

pub struct DatabaseConnection {}

impl DatabaseConnection {
    fn get_database_url() -> String {
        match dotenv::dotenv() {
            Ok(_) => (),
            Err(err) => panic!("proteomic::utility::database_connection::DatabaseConnection::get_database_connection(): Could not load .env-file, reason: {}", err)
        }
        return env::var("PGSQL_URL").expect("proteomic::utility::database_connection::DatabaseConnection::get_database_connection(): Variable 'PGSQL_URL' must be set in .env");
    }

    pub fn get_database_connection() -> postgres::Connection {
        match postgres::Connection::connect(DatabaseConnection::get_database_url().as_str(), postgres::TlsMode::None) {
            Ok(connection) => return connection,
            Err(err) => panic!("proteomic::utility::database_connection::DatabaseConnection::get_database_connection(): Could not connect to database, reason: {}", handle_postgres_error(&err))
        }
    }
}