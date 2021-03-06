use std::fs::File;
use std::io::BufReader;
use std::io::prelude::*;
use std::sync::{Arc, Mutex};

use threadpool::ThreadPool;

use proteomic::models::enzyms::digest_enzym::DigestEnzym;
use proteomic::utility::input_file_digester::file_digester::FileDigester;
use proteomic::utility::database_connection::DatabaseConnection;
use proteomic::models::protein::Protein;
use proteomic::utility::logger::async_queued_logger::AsyncQueuedLogger;
use proteomic::utility::logger::async_performance_logger::AsyncPerformanceLogger;
use proteomic::models::enzyms;




pub struct FastaDigester {
    fasta_file_path: String,
    thread_count: usize,
    max_number_of_missed_cleavages: u8,
    min_peptide_length: usize,
    max_peptide_length: usize,
    message_logger: Arc<AsyncQueuedLogger>,
    unsuccessful_protein_logger: Arc<AsyncQueuedLogger>,
    performance_logger: Arc<Mutex<AsyncPerformanceLogger>>
}

// <E: DigestEnzym + Clone + Send + 'static>
impl FileDigester for FastaDigester {
    fn new(file_path: &str, thread_count: usize, max_number_of_missed_cleavages: u8, min_peptide_length: usize, max_peptide_length: usize) -> FastaDigester {
        return FastaDigester {
            fasta_file_path: file_path.to_owned(),
            thread_count: thread_count,
            max_number_of_missed_cleavages: max_number_of_missed_cleavages,
            min_peptide_length: min_peptide_length,
            max_peptide_length: max_peptide_length,
            message_logger: Arc::new(AsyncQueuedLogger::new("./digest.log")),
            unsuccessful_protein_logger: Arc::new(AsyncQueuedLogger::new("./unsuccessful_proteins.log.fasta")),
            performance_logger: Arc::new(Mutex::new(AsyncPerformanceLogger::new("./digest_performance.csv")))
        }
    }

    fn process_file(&mut self, enzym_name: &str) -> f64 {
        let thread_pool = ThreadPool::new(self.thread_count);

        // open fasta file
        let fasta_file = File::open(&self.fasta_file_path).expect("fasta file not found");
        let fasta_file = BufReader::new(fasta_file);

        let mut header: String = String::new();
        let mut aa_sequence = String::new();

        // This should be dynamicaly updated in the future, depending on the type and number of returned postgresql errors from digestion.
        // When the number of inserted proteins and peptides get high enough (numbers are depending on the capabilities of the database server),
        // the database will return out of memory errors and/or some other errors, because it needs more and more time and memory to process the data.
        // To reduce such errors, the number of peptides per transaction could be reduced.
        let transaction_size = 100;

        // start

        match self.performance_logger.lock() {
            Ok(mut logger) => logger.start_logging(),
            Err(_) => println!("proteomic::utility::input_file_digester::fasta_digester.process_file(): Try to lock poisened mutex for performance logger, performance log not working")
        }

        let start_time: f64 = time::precise_time_s();
        for line in fasta_file.lines() {
            // trim
            let string_line = line.unwrap().as_mut_str().trim().to_owned();
            if !string_line.starts_with(">") {
                aa_sequence.push_str(&string_line);
            } else {
                if header.len() > 0 {
                    let mut protein: Protein = Protein::new(header.as_str(), aa_sequence.as_str());
                    // clone logger pointer for thread
                    let message_logger_ptr = self.message_logger.clone();
                    let unsuccessful_protein_logger_ptr = self.unsuccessful_protein_logger.clone();
                    let mut performance_logger_ptr = self.performance_logger.clone();
                    // clone primitves for thread
                    let max_number_of_missed_cleavages = self.max_number_of_missed_cleavages;
                    let min_peptide_length = self.min_peptide_length;
                    let max_peptide_length = self.max_peptide_length;
                    let enzym_name_clone = enzym_name.to_owned();
                    while thread_pool.queued_count() > 0 {} // prevent flooding the queue with threads, wait that queue is empty before adding new thread
                    thread_pool.execute(move||{
                        let db_conn = DatabaseConnection::get_database_connection();
                        let db_conn_ref = &db_conn;
                        {

                            //let mut enzym = Trypsin::new(db_conn_ref, max_number_of_missed_cleavages, min_peptide_length, max_peptide_length);
                            let mut enzym = enzyms::get(enzym_name_clone.as_str(), db_conn_ref, max_number_of_missed_cleavages, min_peptide_length, max_peptide_length);
                            let summary = enzym.digest(&mut protein, transaction_size);

                            if summary.get_unsolveable_errors_occured() {
                                unsuccessful_protein_logger_ptr.push_back(protein.as_fasta_entry());
                                message_logger_ptr.push_back(summary.get_log());
                            }
                            match performance_logger_ptr.lock() {
                                Ok(logger) => logger.increase_counter_by(
                                    if summary.has_created_protein() { 1 } else { 0 },
                                    summary.get_number_of_created_peptides(),
                                    summary.get_number_of_created_peptide_protein_associations(),
                                    1,
                                    summary.get_number_of_processed_peptides(),
                                    summary.get_number_of_processed_peptide_protein_associations(),
                                ),
                                Err(_) => println!("proteomic::utility::input_file_digester::fasta_digester.process_file(): Try to lock poisened mutex for performance logger, performance log not working")
                            };
                        }
                    });
                    aa_sequence = String::new();
                }
                header = string_line;
            }
        }
        // process last protein
        let mut protein: Protein = Protein::new(header.as_str(), aa_sequence.as_str());
        let db_conn = DatabaseConnection::get_database_connection();
        let db_conn_ref = &db_conn;
        {
            let mut enzym = enzyms::get(enzym_name, db_conn_ref, self.max_number_of_missed_cleavages, self.min_peptide_length, self.max_peptide_length);
            let summary = enzym.digest(&mut protein, transaction_size);
            if summary.get_unsolveable_errors_occured() {
                self.unsuccessful_protein_logger.push_back(protein.as_fasta_entry());
                self.message_logger.push_back(summary.get_log());
            }
            match self.performance_logger.lock() {
                Ok(logger) => logger.increase_counter_by(
                    if summary.has_created_protein() { 1 } else { 0 },
                    summary.get_number_of_created_peptides(),
                    summary.get_number_of_created_peptide_protein_associations(),
                    1,
                    summary.get_number_of_processed_peptides(),
                    summary.get_number_of_processed_peptide_protein_associations(),
                ),
                Err(_) => println!("proteomic::utility::input_file_digester::fasta_digester.process_file(): Try to lock poisened mutex for performance logger, performance log not working")
            };
        }
        // wait for threads
        thread_pool.join();
        let stop_time: f64 = time::precise_time_s();

        return stop_time - start_time;
    }
}

