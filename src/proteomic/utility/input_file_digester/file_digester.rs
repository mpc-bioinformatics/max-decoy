extern crate dotenv;

use std::env;
use std::io::prelude::*;
use std::path::Path;
use std::fs::remove_file;
use std::fs::File;

use proteomic::utility::enzyms::digest_enzym::DigestEnzym;

const START_LINE_FILE_PATH: &str = "./start_line.txt";

pub trait FileDigester<E: DigestEnzym + Clone + Send + 'static> {
    fn new(file_path: &str, thread_count: usize, max_number_of_missed_cleavages: i16, min_peptide_length: usize, max_peptide_length: usize) -> Self;
    fn process_file(&self);
    fn get_database_url() -> String {
        dotenv::dotenv().ok();
        return env::var("PGSQL_URL").expect("Postgresql-Database-URL 'PGSQL_URL' must be set ");
    }

    fn start_line_file_exists() -> bool {
        return Path::new(START_LINE_FILE_PATH).is_file();
    }

    fn get_start_line() -> usize {
        if Self::start_line_file_exists() {
            let mut start_line_file = File::open(START_LINE_FILE_PATH).unwrap();
            let mut contents = String::new();
            start_line_file.read_to_string(&mut contents);
            return match contents.as_str().parse::<usize>() {
                Ok(value) => value,
                Err(_err) => 0
            }
        }
        return 0;
    }

    fn update_start_line_file(line_number: usize) -> bool {
        if Self::start_line_file_exists() {
            let mut file = File::create(START_LINE_FILE_PATH).unwrap();
            return match file.write_all(format!("{}", line_number).as_bytes()) {
                Ok(_value) => true,
                Err(_err) => false
            };
        }
        return false;
    }

    fn remove_start_line_file() -> bool {
        if Self::start_line_file_exists() {
            return match remove_file(START_LINE_FILE_PATH) {
                Ok(_value) => true,
                Err(_err) => false
            };
        }
        return false;
    }
}
