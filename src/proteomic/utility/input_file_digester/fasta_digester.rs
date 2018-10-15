extern crate postgres;
extern crate time;
extern crate threadpool;

use std::env;
use std::fs::File;
use std::io::BufReader;
use std::io::prelude::*;
use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::collections::HashSet;

use self::threadpool::ThreadPool;

use proteomic::utility::enzyms::digest_enzym::DigestEnzym;
use proteomic::utility::input_file_digester::file_digester::FileDigester;
use proteomic::models::persistable::Persistable;
use proteomic::models::protein::Protein;



pub struct FastaDigester<E: DigestEnzym + Clone + Send> {
    fasta_file_path: String,
    enzym: E,
    thread_count: usize
}


impl<E: DigestEnzym + Clone + Send + 'static> FileDigester<E> for FastaDigester<E> {
    fn new(file_path: &str, thread_count: usize, number_of_missed_cleavages: i16, min_peptide_length: usize, max_peptide_length: usize) -> FastaDigester<E> {
        return FastaDigester {
            fasta_file_path: file_path.to_owned(),
            enzym: E::new(number_of_missed_cleavages, min_peptide_length, max_peptide_length),
            thread_count: thread_count
        }
    }

    fn process_file(&self) {

        let thread_pool = ThreadPool::new(self.thread_count);

        let mut transaction_config = postgres::transaction::Config::new();
        transaction_config.isolation_level(postgres::transaction::IsolationLevel::ReadUncommitted);
        transaction_config.read_only(false);
        transaction_config.deferrable(false);


        let fasta_file = File::open(&self.fasta_file_path).expect("fasta file not found");
        let fasta_file = BufReader::new(fasta_file);
        println!("#lines in fasta file = {}...", fasta_file.lines().count());

        let fasta_file = File::open(&self.fasta_file_path).expect("fasta file not found");
        let fasta_file = BufReader::new(fasta_file);

        let mut current_line: usize = 1;
        let start_line: usize = Self::get_start_line();
        println!("start at line {}...", start_line);

        let mut header: String = String::new();
        let mut aa_sequence = String::new();

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
                    let enzym_clone: E = self.enzym.clone();
                    while thread_pool.queued_count() > 0 {} // wait for free resources
                    thread_pool.execute(move||{
                        let database_connection: postgres::Connection = postgres::Connection::connect(Self::get_database_url().as_str(), postgres::TlsMode::None).unwrap();

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
                        overall_peptide_counter_clone.fetch_add(enzym_clone.digest(&database_connection, &mut protein), Ordering::Relaxed);
                    });
                    // update_start_line_file(current_line);
                    aa_sequence = String::new();
                }
                header = string_line;
            }
            current_line += 1;
        }
        // process last protein
        let database_connection: postgres::Connection = postgres::Connection::connect(Self::get_database_url().as_str(), postgres::TlsMode::None).unwrap();
        let mut protein: Protein = Protein::new(header.clone(), aa_sequence);
        protein.save(&database_connection);
        let stop_time: f64 = time::precise_time_s();;
    }
}

