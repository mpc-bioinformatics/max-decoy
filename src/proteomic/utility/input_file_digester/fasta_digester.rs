extern crate postgres;
extern crate time;
extern crate threadpool;

use std::fs::File;
use std::io::BufReader;
use std::io::prelude::*;
use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::collections::HashSet;

use self::threadpool::ThreadPool;

use proteomic::utility::enzyms::digest_enzym::DigestEnzym;
use proteomic::utility::input_file_digester::file_digester::FileDigester;
use proteomic::utility::database_connection::DatabaseConnection;
use proteomic::models::persistable::Persistable;
use proteomic::models::protein::Protein;
use proteomic::utility::logger::async_queued_logger::AsyncQueuedLogger;
use proteomic::utility::logger::simple_logger::SimpleLogger;

const PERFORMANCE_LOG_TIME_STEPS_IN_SECONDS: i64 = 900;


pub struct FastaDigester<E: DigestEnzym + Clone + Send> {
    fasta_file_path: String,
    enzym: E,
    thread_count: usize,
    start_line: usize,
    message_logger: Arc<AsyncQueuedLogger>,
    unsuccessful_protein_logger: Arc<AsyncQueuedLogger>,
    performance_logger: SimpleLogger
}


impl<E: DigestEnzym + Clone + Send + 'static> FileDigester<E> for FastaDigester<E> {
    fn new(file_path: &str, thread_count: usize, number_of_missed_cleavages: i16, min_peptide_length: usize, max_peptide_length: usize) -> FastaDigester<E> {
        return FastaDigester {
            fasta_file_path: file_path.to_owned(),
            enzym: E::new(number_of_missed_cleavages, min_peptide_length, max_peptide_length),
            thread_count: thread_count,
            start_line: Self::get_start_line(),
            message_logger: Arc::new(AsyncQueuedLogger::new("./digest.log")),
            unsuccessful_protein_logger: Arc::new(AsyncQueuedLogger::new("./unsuccessful_proteins.log.fasta")),
            performance_logger: SimpleLogger::new("./digest_performance.csv")
        }
    }

    fn process_file(&mut self) -> (usize, usize, f64) {
        let thread_pool = ThreadPool::new(self.thread_count); // TODO: do not assign this if self.thread_count < 2

        // database coonnection for main thread
        let database_connection: postgres::Connection = DatabaseConnection::get_database_connection();

        // open fasta file for counting lines
        // let fasta_file = File::open(&self.fasta_file_path).expect("fasta file not found");
        // let fasta_file = BufReader::new(fasta_file);
        // println!("#lines in fasta file = {}...", fasta_file.lines().count());

        // open fasta file
        let fasta_file = File::open(&self.fasta_file_path).expect("fasta file not found");
        let fasta_file = BufReader::new(fasta_file);

        // line counter
        let mut current_line: usize = 1;
        println!("start at line {}...", self.start_line);

        let mut header: String = String::new();
        let mut aa_sequence = String::new();

        // thread safe counter
        let commited_protein_counter = Arc::new(AtomicUsize::new(0));
        let commited_peptide_counter = Arc::new(AtomicUsize::new(0));
        let commited_peptide_protein_association_counter = Arc::new(AtomicUsize::new(0));

        // performance measurement
        self.performance_logger.write(&String::from("\"seconds\",\"proteins\",\"peptides\",\"protein/peptide-association\"\n"));
        let start_at_sec = time::now().to_timespec().sec;
        let mut next_performance_log_sec = start_at_sec + PERFORMANCE_LOG_TIME_STEPS_IN_SECONDS;

        let start_time: f64 = time::precise_time_s();
        for line in fasta_file.lines() {
            if current_line < self.start_line {
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
                    if self.thread_count > 1 {
                        // clones for thread
                        let mut commited_protein_counter_clone = commited_protein_counter.clone();
                        let mut commited_peptide_counter_clone = commited_peptide_counter.clone();
                        let mut commited_peptide_protein_association_counter_clone = commited_peptide_protein_association_counter.clone();
                        let enzym_clone: E = self.enzym.clone();
                        let message_logger_clone = self.message_logger.clone();
                        let unsuccessful_protein_logger_clone = self.unsuccessful_protein_logger.clone();
                        while thread_pool.queued_count() > 0 {} // prevent flooding the queue with threads, wait that queue is empty before adding new thread
                        thread_pool.execute(move||{
                            let database_connection: postgres::Connection = DatabaseConnection::get_database_connection();
                            match protein.save(&database_connection) {
                                Ok(_) => {
                                    commited_protein_counter_clone.fetch_add(1, Ordering::Relaxed);
                                    match enzym_clone.digest(&database_connection, &mut protein){
                                        Ok(digest_ok) => {
                                            if digest_ok.get_log_message().len() > 0 {
                                                message_logger_clone.as_ref().push_back(digest_ok.get_log_message().to_owned());
                                            }
                                            commited_peptide_counter_clone.fetch_add(digest_ok.get_commited_peptides(), Ordering::Relaxed);
                                            commited_peptide_protein_association_counter_clone.fetch_add(digest_ok.get_commited_peptide_protein_associations(), Ordering::Relaxed);
                                            ()
                                        },
                                        Err(digest_error) => {
                                            unsuccessful_protein_logger_clone.as_ref().push_back(protein.as_fasta_entry());
                                            message_logger_clone.as_ref().push_back(digest_error.get_message().to_owned());
                                            ()
                                        }
                                    }
                                },
                                Err(err) =>{
                                    message_logger_clone.as_ref().push_back(format!("ERROR   => PROTEIN [{}]: Could not be commited, error\n\t{}\n", protein.get_aa_sequence(), err));
                                    unsuccessful_protein_logger_clone.as_ref().push_back(protein.as_fasta_entry());
                                }
                            }
                        });
                    } else {
                        match protein.save(&database_connection) {
                            Ok(_) => {
                                commited_protein_counter.fetch_add(1, Ordering::Relaxed);
                                match self.enzym.digest(&database_connection, &mut protein) {
                                    Ok(digest_ok) => {
                                        if digest_ok.get_log_message().len() > 0 {
                                            self.message_logger.as_ref().push_back(digest_ok.get_log_message().to_owned());
                                        }
                                        commited_peptide_counter.fetch_add(digest_ok.get_commited_peptides(), Ordering::Relaxed);
                                        commited_peptide_protein_association_counter.fetch_add(digest_ok.get_commited_peptide_protein_associations(), Ordering::Relaxed);
                                        ()
                                    },
                                    Err(digest_error) => {
                                        self.unsuccessful_protein_logger.push_back(protein.as_fasta_entry());
                                        self.message_logger.push_back(digest_error.get_message().to_owned());
                                        ()
                                    }
                                }
                            },
                            Err(err) =>{
                                self.message_logger.as_ref().push_back(format!("ERROR   => PROTEIN [{}]: Could not be commited, error\n\t{}\n", protein.get_aa_sequence(), err));
                                self.unsuccessful_protein_logger.as_ref().push_back(protein.as_fasta_entry());
                            }
                        }
                    }
                    // update_start_line_file(current_line);
                    aa_sequence = String::new();
                    let now_in_seconds = time::now().to_timespec().sec;
                    if now_in_seconds >= next_performance_log_sec {
                        self.performance_logger.write(
                            &format!(
                                "{},{},{},{}\n",
                                now_in_seconds - start_at_sec,
                                commited_protein_counter.load(Ordering::Relaxed),
                                commited_peptide_counter.load(Ordering::Relaxed),
                                commited_peptide_protein_association_counter.load(Ordering::Relaxed),
                            )
                        );
                        next_performance_log_sec = now_in_seconds + PERFORMANCE_LOG_TIME_STEPS_IN_SECONDS;
                    }
                }
                header = string_line;
            }
            current_line += 1;
        }
        // process last protein
        let mut protein: Protein = Protein::new(header.clone(), aa_sequence);
        match protein.save(&database_connection) {
            Ok(_) => {
                commited_protein_counter.fetch_add(1, Ordering::Relaxed);
                match self.enzym.digest(&database_connection, &mut protein) {
                    Ok(digest_ok) => {
                        if digest_ok.get_log_message().len() > 0 {
                            self.message_logger.as_ref().push_back(digest_ok.get_log_message().to_owned());
                        }
                        commited_peptide_counter.fetch_add(digest_ok.get_commited_peptides(), Ordering::Relaxed);
                        commited_peptide_protein_association_counter.fetch_add(digest_ok.get_commited_peptide_protein_associations(), Ordering::Relaxed);
                        ()
                    },
                    Err(digest_error) => {
                        self.unsuccessful_protein_logger.push_back(protein.as_fasta_entry());
                        self.message_logger.push_back(digest_error.get_message().to_owned());
                        ()
                    }
                }
            },
            Err(err) =>{
                self.message_logger.as_ref().push_back(format!("ERROR   => PROTEIN [{}]: Could not be commited, error\n\t{}\n", protein.get_aa_sequence(), err));
                self.unsuccessful_protein_logger.as_ref().push_back(protein.as_fasta_entry());
            }
        }
        // wait for threads
        thread_pool.join();
        let stop_time: f64 = time::precise_time_s();
        let now_in_seconds = time::now().to_timespec().sec;
        self.performance_logger.write(
            &format!(
                "{},{},{},{}\n",
                now_in_seconds - start_at_sec,
                commited_protein_counter.load(Ordering::Relaxed),
                commited_peptide_counter.load(Ordering::Relaxed),
                commited_peptide_protein_association_counter.load(Ordering::Relaxed),
            )
        );
        return (commited_protein_counter.load(Ordering::Relaxed), commited_peptide_counter.load(Ordering::Relaxed), stop_time - start_time);
    }

    fn process_file_but_count_only(&self) -> (usize, usize) {
        let fasta_file = File::open(&self.fasta_file_path).expect("fasta file not found");
        let fasta_file = BufReader::new(fasta_file);

        // line counter
        let mut current_line: usize = 1;

        // vars
        let mut protein_counter: usize = 0;
        let mut header: String = String::new();
        let mut aa_sequence = String::new();
        let mut set_of_aa_sequences: HashSet<String> = HashSet::new();

        for line in fasta_file.lines() {
            if current_line < self.start_line {
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
                    protein_counter += 1;
                    self.enzym.digest_with_hash_set(&mut protein, &mut set_of_aa_sequences);
                    aa_sequence = String::new();
                }
                header = string_line;
            }
        }
        let mut protein: Protein = Protein::new(header.clone(), aa_sequence);
        protein_counter += 1;
        self.enzym.digest_with_hash_set(&mut protein, &mut set_of_aa_sequences);
        return (protein_counter, set_of_aa_sequences.len());
    }
}

