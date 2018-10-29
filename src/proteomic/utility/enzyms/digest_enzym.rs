extern crate onig;
extern crate postgres;

use std::{thread, time};
use std::collections::HashSet;

use proteomic::models::persistable::Persistable;
use proteomic::models::protein::Protein;
use proteomic::models::peptide::Peptide;
use proteomic::models::peptide_protein_association::PeptideProteinAssociation;
use proteomic::utility::enzyms::results::digest_ok::DigestOk;
use proteomic::utility::enzyms::results::digest_error::DigestError;

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

            // database transaction and statements
            let transaction = database_connection.transaction().unwrap();

            let peptide_insert_statement = transaction.prepare_cached(Peptide::get_insert_query()).unwrap();
            let peptide_select_by_unique_identifier_statement = transaction.prepare_cached(Peptide::get_select_primary_key_by_unique_identifier_query()).unwrap();

            let peptide_protein_association_insert_statement = transaction.prepare_cached(PeptideProteinAssociation::get_insert_query()).unwrap();
            let peptide_protein_association_select_by_unique_identifier_statement = transaction.prepare_cached(PeptideProteinAssociation::get_select_primary_key_by_unique_identifier_query()).unwrap();
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
                            processed_peptide_counter += 1;
                            error_message = String::new();
                            let mut peptide = Peptide::new(new_peptide_aa_sequence.clone(), self.get_shortcut().to_owned(), number_of_missed_cleavages);
                            match peptide.exec_insert_statement(&peptide_insert_statement) {
                                Ok(_) => commited_peptide_counter += 1,
                                Err(insert_err) => match insert_err.as_str() {
                                    // no return means the peptide already exists, this
                                    "NORET" => {
                                        match peptide.exec_select_primary_key_by_unique_identifier_statement(&peptide_select_by_unique_identifier_statement) {
                                            Ok(_) => (),
                                            // ok, error at this point means the peptides could not be insertet nor selected absolutely an error
                                            Err(select_err) => {
                                                error_counter += 1;
                                                error_message.push_str(format!("\n\tinsert: {}\n\tselect: {}", insert_err, select_err).as_str());
                                                break 'peptide_loop;
                                            }
                                        }
                                    },
                                    // anything else than NORET is bad
                                    _ => {
                                        error_message.push_str(format!("\n\tinsert: {}", insert_err).as_str());
                                        error_counter += 1;
                                        break 'peptide_loop;
                                    }
                                }
                            }
                            if peptide.is_persisted() {
                                processed_peptide_protein_association_counter += 1;
                                let mut association = PeptideProteinAssociation::new(&peptide, &protein);
                                match association.exec_insert_statement(&peptide_protein_association_select_by_unique_identifier_statement) {
                                    Ok(_) => processed_peptide_protein_association_counter += 1,
                                    Err(_insert_err) => {
                                        match association.exec_select_primary_key_by_unique_identifier_statement(&peptide_protein_association_insert_statement){
                                            Ok(_) => commited_peptide_protein_association_counter += 1,
                                            Err(_select_err) => () //error_message.push_str(format!("select: {}", insert_err).to_str())
                                        }
                                    }
                                }
                            }
                        }
                    } else {
                        break 'missed_cleavage_loop;
                    }
                }
                // peptide_position += peptides_without_missed_cleavages[peptide_idx].len();
            }
            //println!("THREAD [{}]: error_counter => {}", protein.get_accession(), error_counter);
            if !error_message.is_empty() {
                match error_counter {
                    // some errors occure, rollback, wait and try again later
                    1 | 2 => {
                        log.push_str(
                            format!(
                                "WARNING => THREAD [{}]: {}. error occured. Do a rollback and try again in {} seconds. Error(s):{}\n",
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