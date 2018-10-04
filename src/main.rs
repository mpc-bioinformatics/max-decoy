extern crate time;
extern crate postgres;
extern crate dotenv;
extern crate num_cpus;
extern crate threadpool;

use std::env;
use std::fs::File;
use std::io::BufReader;
use std::io::prelude::*;
use std::path::Path;
use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering};

use threadpool::ThreadPool;

mod proteomic;
use proteomic::utility::enzym::{DigestEnzym, Trypsin};
use proteomic::models::persistable::Persistable;
use proteomic::models::protein::Protein;

mod tests;

const START_LINE_FILE_PATH: &str = "./start_line.txt";

fn get_database_url() -> String {
    dotenv::dotenv().ok();
    return env::var("PGSQL_URL").expect("Postgresql-Database-URL 'PGSQL_URL' must be set ");
}

fn start_line_file_exists() -> bool {
    return Path::new(START_LINE_FILE_PATH).is_file();
}

fn get_start_line() -> usize {
    if start_line_file_exists() {
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
    if start_line_file_exists() {
        let mut file = File::create(START_LINE_FILE_PATH).unwrap();
        return match file.write_all(format!("{}", line_number).as_bytes()) {
            Ok(_value) => true,
            Err(_err) => false
        };
    }
    return false;
}

fn remove_start_line_file() -> bool {
    if start_line_file_exists() {
        return match std::fs::remove_file(START_LINE_FILE_PATH) {
            Ok(_value) => true,
            Err(_err) => false
        };
    }
    return false;
}

fn main() {
    let args: Vec<String> = env::args().collect();

    let filename = &args[1];
    println!("use fasta file {}...", filename);

    let cpus: usize = num_cpus::get() - 1;
    let thread_pool = ThreadPool::new(cpus);
    println!("use {} threads...", cpus + 1);


    let fasta_file = File::open(filename).expect("fasta file not found");
    let fasta_file = BufReader::new(fasta_file);
    println!("#lines in fasta file = {}...", fasta_file.lines().count());

    let fasta_file = File::open(filename).expect("fasta file not found");
    let fasta_file = BufReader::new(fasta_file);

    let mut current_line: usize = 1;
    let start_line: usize = get_start_line();
    println!("start at line {}...", start_line);

    let mut header: String = String::new();
    let mut aa_sequence = String::new();
    let trypsin: Trypsin = Trypsin::new(2, 6, 50);


    let overall_protein_counter = Arc::new(AtomicUsize::new(0));
    let overall_peptide_counter = Arc::new(AtomicUsize::new(0));

    let start_time: f64 = time::precise_time_s();
    for line in fasta_file.lines() {
        if current_line < start_line {
            current_line += 1;
            continue;
        }
        // trim
        let string_line = line.unwrap().as_mut_str().trim().to_owned();
        if !string_line.starts_with(">") {
            aa_sequence.push_str(&string_line);
        } else {
            if header.len() > 0 {
                let mut overall_protein_counter_clone = overall_protein_counter.clone();
                let mut overall_peptide_counter_clone = overall_peptide_counter.clone();
                let trypsin_clone: Trypsin = trypsin.clone();
                while thread_pool.queued_count() > 0 {} // wait for free resources
                thread_pool.execute(move||{
                    let database_connection: postgres::Connection = postgres::Connection::connect(get_database_url().as_str(), postgres::TlsMode::None).unwrap();
                    let mut protein: Protein = Protein::new(header.clone(), aa_sequence);
                    if !protein.save(&database_connection) {
                        println!("protein not saved:\n\t{}\n", protein.to_string());
                    }
                    overall_protein_counter_clone.fetch_add(1, Ordering::Relaxed);
                    overall_peptide_counter_clone.fetch_add(trypsin_clone.digest(&database_connection, &mut protein), Ordering::Relaxed);
                });
                // update_start_line_file(current_line);
                aa_sequence = String::new();
            }
            header = string_line;
        }
        current_line += 1;
    }
    // process last protein
    let database_connection: postgres::Connection = postgres::Connection::connect(get_database_url().as_str(), postgres::TlsMode::None).unwrap();
    let mut protein: Protein = Protein::new(header.clone(), aa_sequence);
    protein.save(&database_connection);
    overall_protein_counter.fetch_add(1, Ordering::Relaxed);
    overall_peptide_counter.fetch_add(trypsin.digest(&database_connection, &mut protein), Ordering::Relaxed);
    // remove_start_line_file();
    let stop_time: f64 = time::precise_time_s();
    println!("Proteins processed: {}", overall_protein_counter.load(Ordering::Relaxed));
    println!("  Peptides created: {}", overall_peptide_counter.load(Ordering::Relaxed));
    println!("Need {} s", (stop_time - start_time));
}

#[cfg(test)]
fn test() {
    tests::run();
}