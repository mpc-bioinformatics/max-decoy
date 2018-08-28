extern crate postgres;
extern crate dotenv;

pub mod protein;
pub mod peptide;
pub mod peptide_position;
pub mod collection;
pub mod persistable;

use std::env;


pub fn get_db_connection() -> postgres::Connection {
    dotenv::dotenv().ok();
    let database_url = env::var("PGSQL_URL").expect("Postgresql-Database-URL 'PGSQL_URL' must be set ");
    // let tls_env = env::var("PGSQL_TLS_MODE").expect("Postgresql TLS-mode 'PGSQL_TLS_MODE' must be set");
    // let tls_mode: postgres::TlsMode = match &tls_env.as_ref() {
    //     "NONE" => postgres::TlsMode::None,
    //     // "PREFER" => "postgres::TlsMode::Prefer",
    //     // "REQUIRE" => "postgres::TlsMode::Require",
    //     // _ => "TlsMode::None",
    // };
    return postgres::Connection::connect(database_url, postgres::TlsMode::None).unwrap();
}