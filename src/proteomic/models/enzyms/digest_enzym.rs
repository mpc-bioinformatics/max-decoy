extern crate onig;
extern crate postgres;

use std::{thread, time};
use std::collections::HashSet;

use proteomic::models::persistable::{Persistable, QueryError, QueryOk};
use proteomic::models::protein::Protein;
use proteomic::models::peptide::Peptide;
use proteomic::models::peptide_protein_association::PeptideProteinAssociation;
use proteomic::models::enzyms::results::digest_ok::DigestOk;
use proteomic::models::enzyms::results::digest_error::DigestError;

const DIGEST_WAIT_DURATION_FOR_ERRORS: time::Duration = time::Duration::from_secs(20);

pub trait DigestEnzym {
    fn new(max_number_of_missed_cleavages: i16, min_peptide_length: usize, max_peptide_length: usize) -> Self;
    fn get_name(&self) -> &str;
    fn get_shortcut(&self) -> &str;
    fn get_max_number_of_missed_cleavages(&self) -> i16;
    fn get_digest_regex(&self) -> &onig::Regex;
    fn get_digest_replace(&self) -> &'static str;
    fn get_min_peptide_length(&self) -> usize;
    fn get_max_peptide_length(&self) -> usize;

    fn digest(&self, database_connection: &postgres::Connection, protein: &mut Protein) -> Result<DigestOk, DigestError>  {
        let mut error_counter: u8 = 0;
        let mut commited_peptide_counter: usize = 0;
        let mut commited_peptide_protein_association_counter: usize = 0;
        let mut processed_peptide_counter: usize = 0;
        let mut processed_peptide_protein_association_counter: usize = 0;
        let mut log: String = String::new();

        /*
         * clone aa_squence and pass it as mutable into replace_all
         * replace every digist_regex-match with with digist_replace (in caseof Trypsin it means add a whitespace between K or T and not P)
         * make the result mutable and split it on whitespaces
         * collect the results as String-vector
         */
        let peptides_without_missed_cleavages: Vec<String> = self.get_digest_regex().split(protein.get_aa_sequence().clone().as_mut_str()).map(|peptide| peptide.to_owned()).collect::<Vec<String>>();

        'tries_loop: loop {
            processed_peptide_counter = 0;
            processed_peptide_protein_association_counter = 0;
            commited_peptide_counter = 0;                        // reset peptides counter for this try
            commited_peptide_protein_association_counter = 0;    // reset association  counter for this try

            let mut error_message: String = String::new();
            let mut error_occured = false;

            // database transaction and statements
            let transaction = database_connection.transaction().unwrap();

            let peptide_create_statement = transaction.prepare_cached(Peptide::create_query()).unwrap();
            let peptide_exists_statement = transaction.prepare_cached(Peptide::exists_query()).unwrap();

            let pp_association_create_statement = transaction.prepare_cached(PeptideProteinAssociation::create_query()).unwrap();
            let pp_association_exists_statement = transaction.prepare_cached(PeptideProteinAssociation::exists_query()).unwrap();
            ////

            // let mut peptide_position: usize = 0;
            // calculate peptides for missed_cleavages 1 to n + 1 (+1 because explicit boundary)
            'peptide_loop: for peptide_idx in 0..peptides_without_missed_cleavages.len() {
                let mut new_peptide_aa_sequence: String = String::new();
                'missed_cleavage_loop: for number_of_missed_cleavages in 0..(self.get_max_number_of_missed_cleavages() + 1) {
                    let temp_idx: usize = peptide_idx + number_of_missed_cleavages as usize;
                    if temp_idx < peptides_without_missed_cleavages.len() {
                        new_peptide_aa_sequence.push_str(peptides_without_missed_cleavages.get(temp_idx).unwrap());
                        if self.is_aa_sequence_in_range(&new_peptide_aa_sequence) {
                            let mut peptide = Peptide::new(new_peptide_aa_sequence.clone(), self.get_shortcut().to_owned(), number_of_missed_cleavages);
                            // reset error variables
                            error_message = format!("Peptide [{}]", peptide.get_aa_sequence());
                            error_occured = false;
                            match peptide.prepared_create(&peptide_create_statement, &peptide_exists_statement) {
                                Ok(query_ok) => match query_ok {
                                    QueryOk::Created => commited_peptide_counter += 1,
                                    QueryOk::AlreadyExists => (),
                                    _ => panic!("Panic [proteomic::models::enzyms::digest_enzym::DigestEnzym::digest()]: In fact no other QueryOk than QueryOk::Created and QueryOk::AlreadyExists used in Peptide.prepared_create(), this panic should never be reached")
                                },
                                Err(err) => {
                                    error_occured = true;
                                    error_message.push_str(format!("\n\tPeptide.prepared_create(): {}", err).as_str());
                                    break 'peptide_loop;
                                }
                            }
                            processed_peptide_counter += 1;
                            if peptide.is_persisted() {
                                let mut association = PeptideProteinAssociation::new(&peptide, &protein);
                                match association.prepared_create(&pp_association_create_statement, &pp_association_exists_statement) {
                                    Ok(query_ok) => match query_ok {
                                        QueryOk::Created => commited_peptide_protein_association_counter += 1,
                                        QueryOk::AlreadyExists => (),
                                        _ => panic!("Panic [proteomic::models::enzyms::digest_enzym::DigestEnzym::digest()]: In fact no other QueryOk than QueryOk::Created and QueryOk::AlreadyExists used in PeptideProteinAssociation.prepared_create(), this panic should never be reached")
                                    },
                                    Err(err) => {
                                        error_occured = true;
                                        error_message.push_str(format!("\n\tPeptideProteinAssociation.prepared_create(): {}",err).as_str());
                                        break 'peptide_loop;
                                    }                                     
                                }
                                processed_peptide_protein_association_counter += 1;
                            }
                        }
                    } else {
                        break 'missed_cleavage_loop;
                    }
                }
            }
            if error_occured {
                error_counter += 1;
                match error_counter {
                    // some errors occure, rollback, wait and try again later
                    1 | 2 => {
                        log.push_str(
                            format!(
                                "WARNING => THREAD [{}]: {}. error occured. Do a rollback and try again in {} seconds. Error(s): {}\n",
                                protein.get_accession(),
                                error_counter,
                                DIGEST_WAIT_DURATION_FOR_ERRORS.as_secs(),
                                error_message
                            ).as_str()
                        );
                        transaction.set_rollback();
                        transaction.finish();
                        thread::sleep(DIGEST_WAIT_DURATION_FOR_ERRORS);
                        continue 'tries_loop;
                    },
                    // no more try, rollback, report this as error and set peptides counter to 0
                    _ => {
                        log.push_str(
                            format!(
                                "ERROR   => THREAD [{}]: {}. error occured. Do a rollback, do no further tries. Last occured error(s):{}\n",
                                protein.get_accession(),
                                error_counter,
                                error_message
                            ).as_str()
                        );
                        transaction.set_rollback();
                        transaction.finish();
                        return Err(
                            DigestError::new(
                                log.as_str()
                            )
                        )
                    }

                }
            } else {
                // no errors, commit
                transaction.set_commit();
                transaction.finish();
                if error_counter > 0 {
                    log.push_str(
                        format!(
                            "INFO    => THREAD [{}]: Commited after {} errors occured.\n",
                            protein.get_accession(),
                            error_counter
                        ).as_str()
                    );
                }
                return Ok(
                    DigestOk::new(
                        processed_peptide_counter,
                        commited_peptide_counter,
                        processed_peptide_protein_association_counter,
                        commited_peptide_protein_association_counter,
                        log.as_str()
                    )
                )
            }
        }
    }

    fn digest_with_hash_set(&self, protein: &mut Protein, aa_sequences: &mut HashSet<String>) {
        /*
         * clone aa_squence and pass it as mutable into replace_all
         * replace every digist_regex-match with with digist_replace (in caseof Trypsin it means add a whitespace between K or T and not P)
         * make the result mutable and split it on whitespaces
         * collect the results as String-vector
         */
        let peptides_without_missed_cleavages: Vec<String> = self.get_digest_regex().split(protein.get_aa_sequence().clone().as_mut_str()).map(|peptide| peptide.to_owned()).collect::<Vec<String>>();
        // let mut peptide_position: usize = 0;
        // calculate peptides for missed_cleavages 1 to n + 1 (+1 because explicit boundary)
        'outer: for peptide_idx in 0..peptides_without_missed_cleavages.len() {
            let mut new_peptide_aa_sequence: String = String::new();
            'inner: for number_of_missed_cleavages in 0..(self.get_max_number_of_missed_cleavages() + 1) {
                let temp_idx: usize = peptide_idx + number_of_missed_cleavages as usize;
                if temp_idx < peptides_without_missed_cleavages.len() {
                    new_peptide_aa_sequence.push_str(peptides_without_missed_cleavages.get(temp_idx).unwrap());
                    if self.is_aa_sequence_in_range(&new_peptide_aa_sequence) {
                        aa_sequences.insert(Peptide::gerneralize_aa_sequence(&new_peptide_aa_sequence));
                    }
                } else {
                    break;
                }
            }
        }
    }


    fn is_aa_sequence_in_range(&self, aa_sequence: &String) -> bool {
        return self.get_min_peptide_length() <= aa_sequence.len() && aa_sequence.len() <= self.get_max_peptide_length();
    }
}