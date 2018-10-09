extern crate time;
extern crate postgres;
extern crate dotenv;
extern crate num_cpus;
extern crate threadpool;

use std::env;
use std::fs::File;
use std::io::BufReader;
use std::io::prelude::*;
use std::io::LineWriter;
use std::path::Path;
use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::collections::HashSet;

use threadpool::ThreadPool;

mod proteomic;
use proteomic::utility::enzym::{DigestEnzym, Trypsin};
use proteomic::models::persistable::Persistable;
use proteomic::models::protein::Protein;
use proteomic::models::peptide::Peptide;

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

    // thread safe counter
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
                    match protein.save(&database_connection) {
                        Ok(_) => (),
                        Err(err) => {
                            println!(
                                "ERROR [INSERT & SELECT PROTEIN]:\n\tprotein: {}\n\terror: {}",
                                protein.get_accession(),
                                err
                            );
                        }
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
    let stop_time: f64 = time::precise_time_s();;



    ////////////////////////////////////////////////////////////////////////////
    //// calculate desired values
    // file reader
    let fasta_file = File::open(filename).expect("fasta file not found");
    let fasta_file = BufReader::new(fasta_file);

    // vars
    let mut desired_protein_counter: usize = 0;
    let mut header: String = String::new();
    let mut aa_sequence = String::new();
    let mut desired_aa_sequences: HashSet<String> = HashSet::new();

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
                let mut protein: Protein = Protein::new(header.clone(), aa_sequence);
                desired_protein_counter += 1;
                trypsin.digest_with_hash_set(&mut protein, &mut desired_aa_sequences);
                aa_sequence = String::new();
            }
            header = string_line;
        }
        current_line += 1;
    }
    // process last protein
    let mut protein: Protein = Protein::new(header.clone(), aa_sequence);
    desired_protein_counter += 1;
    trypsin.digest_with_hash_set(&mut protein, &mut desired_aa_sequences);
    ///////////////////////////////////////////////////////////////////////////////////
    //// report
    let actual_number_of_proteins = overall_protein_counter.load(Ordering::Relaxed);
    let actual_number_of_peptides_created_in_threads = overall_peptide_counter.load(Ordering::Relaxed);
    let proteins_in_database = match database_connection.query("SELECT cast(count(id) AS BIGINT) FROM proteins", &[]) {
        Ok(ref rows) if rows.len() > 0 => rows.get(0).get::<usize, i64>(0),
        _ => -1
    };
    let peptides_in_database = match database_connection.query("SELECT cast(count(id) AS BIGINT) FROM peptides", &[]) {
        Ok(ref rows) if rows.len() > 0 => rows.get(0).get::<usize, i64>(0),
        _ => -1
    };


    println!("\nReport:");
    println!("\tseconds needed for calculating and committing: {}", (stop_time - start_time));
    println!(
        "{:<20}{:<20}{:<20}{:<20}{:<20}{:<20}{:<20}{:<20}",
        "type",
        "desired",
        "actual",
        "desired/actual",
        "desired/actual %",
        "in db",
        "desired/db",
        "desired/db %"
    );

    println!(
        "{:<20}{:<20}{:<20}{:<20}{:<20}{:<20}{:<20}{:<20}",
        "proteins",
        desired_protein_counter,
        actual_number_of_proteins,
        desired_protein_counter as i64 - actual_number_of_proteins as i64,
        100.0 / desired_protein_counter as f32 * actual_number_of_proteins as f32,
        proteins_in_database,
        desired_protein_counter as i64 - proteins_in_database,
        100.0 / desired_protein_counter as f32 * proteins_in_database as f32,
    );

    println!(
        "{:<20}{:<20}{:<20}{:<20}{:<20}{:<20}{:<20.4}{:<20}\n\n",
        "peptides",
        desired_aa_sequences.len(),
        actual_number_of_peptides_created_in_threads,
        desired_aa_sequences.len() as i64 - actual_number_of_peptides_created_in_threads as i64,
        100.0 / desired_aa_sequences.len() as f32 * actual_number_of_peptides_created_in_threads as f32,
        peptides_in_database,
        desired_aa_sequences.len() as i64 - peptides_in_database,
        100.0 / desired_aa_sequences.len() as f32 * peptides_in_database as f32,
    );

}

#[cfg(test)]
fn test() {
    tests::run();
}