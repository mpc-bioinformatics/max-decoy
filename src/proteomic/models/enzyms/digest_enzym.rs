use std::{thread, time};
use std::collections::HashSet;

use proteomic::models::persistable::{handle_postgres_error, Persistable, QueryOk};
use proteomic::models::protein::Protein;
use proteomic::models::peptide::Peptide;
use proteomic::models::peptide_protein_association::PeptideProteinAssociation;
use proteomic::models::amino_acids::amino_acid::AminoAcid;
use proteomic::models::enzyms::digest_summary::DigestSummary;

const DIGEST_WAIT_DURATION_FOR_ERRORS: time::Duration = time::Duration::from_secs(5);

pub trait DigestEnzym<'e> {
    fn new(database_connection: &'e  postgres::Connection, max_number_of_missed_cleavages: i16, min_peptide_length: usize, max_peptide_length: usize) -> Self;
    fn get_name(&self) -> &str;
    fn get_shortcut(&self) -> &str;
    fn get_max_number_of_missed_cleavages(&self) -> i16;
    fn get_digest_regex(&self) -> &onig::Regex;
    fn get_digest_replace(&self) -> &'static str;
    fn get_min_peptide_length(&self) -> usize;
    fn get_max_peptide_length(&self) -> usize;
    fn get_database_connection(&self) -> &postgres::Connection;
    fn get_peptide_create_statement(&self) -> &postgres::stmt::Statement;
    fn get_peptide_exists_statement(&self) -> &postgres::stmt::Statement;
    fn get_pp_association_create_statement(&self) -> &postgres::stmt::Statement;
    fn get_pp_association_exists_statement(&self) -> &postgres::stmt::Statement;



    // -> Result<DigestOk, DigestError>
    fn digest(&mut self, protein: &mut Protein) -> DigestSummary  {
        let mut summary = DigestSummary::new();

        // make sure protein is persisted (has id)
        if !protein.is_persisted() {
            match protein.create(self.get_database_connection()) {
                Ok(query_ok) => match query_ok {
                    QueryOk::Created => summary.set_created_protein(),
                    QueryOk::AlreadyExists => (),
                    _ => panic!("proteomic::models::enzyms::digest_enzym::DigestEnzym::digest(): In fact no other QueryOk than QueryOk::Created and QueryOk::AlreadyExists used in Protein.create(), this panic should never be reached")
                },
                // early return, because without persisted Protein we could not create associations
                Err(query_err) => {
                    summary.set_unsolveable_errors_occured();
                    summary.log_push(format!("proteomic::models::enzyms::digest_enzym::DigestEnzym.digest(): Error occured at protein.create(): {}", query_err).as_str());
                    return summary;
                }
            }
        }

        /*
         * clone aa_squence and pass it as mutable into replace_all
         * replace every digist_regex-match with with digist_replace (in caseof Trypsin it means add a whitespace between K or T not followed by P)
         * make the result mutable and split it on whitespaces
         * collect the results as String-vector
         */
        let peptides_without_missed_cleavages: Vec<String> = self.get_digest_regex().split(protein.get_aa_sequence()).map(|peptide| peptide.to_owned()).collect::<Vec<String>>();
        // let mut peptide_position: usize = 0;
        // calculate peptides for missed_cleavages 1 to n + 1 (+1 because explicit boundary)
        for peptide_idx in 0..peptides_without_missed_cleavages.len() {
            // create empty sequence
            let mut new_peptide_aa_sequence: String = String::new();
            'missed_cleavage_loop: for number_of_missed_cleavages in 0..(self.get_max_number_of_missed_cleavages() + 1) {
                let temp_idx: usize = peptide_idx + number_of_missed_cleavages as usize;
                if temp_idx < peptides_without_missed_cleavages.len() {
                    new_peptide_aa_sequence.push_str(peptides_without_missed_cleavages.get(temp_idx).unwrap());
                    if self.is_aa_sequence_in_range(&new_peptide_aa_sequence) {
                        let mut local_log: Vec<String> = Vec::new();
                        let mut peptide = Peptide::new(new_peptide_aa_sequence.as_str(), number_of_missed_cleavages);
                        for try in 1..=3 {
                            match self.create_peptide_and_association(&mut summary, protein, &mut peptide) {
                                Ok(_) => (),
                                Err(err_message) => match try {
                                    1 => {
                                        local_log.push(format!("proteomic::models::enzyms::digest_enzym::DigestEnzym::digest(): transaction for peptide {} failed on 1. try. reason:\n\t{}", peptide.get_aa_sequence(), err_message));
                                        thread::sleep(DIGEST_WAIT_DURATION_FOR_ERRORS);
                                    }
                                    2 =>{
                                        local_log.push(format!("proteomic::models::enzyms::digest_enzym::DigestEnzym::digest(): transaction for peptide {} failed on 2. try. reason:\n\t{}", peptide.get_aa_sequence(), err_message));
                                        thread::sleep(DIGEST_WAIT_DURATION_FOR_ERRORS);
                                    }
                                    3 => {
                                        summary.set_unsolveable_errors_occured();
                                        local_log.push(format!("proteomic::models::enzyms::digest_enzym::DigestEnzym::digest(): transaction for peptide {} failed on 3. try. reason:\n\t{}", peptide.get_aa_sequence(), err_message));
                                        summary.log_push(local_log.join("\n").as_str());
                                    }
                                    _ => ()
                                }
                            }
                        }
                    }
                } else {
                    break 'missed_cleavage_loop;
                }
            }
        }
        if !summary.get_unsolveable_errors_occured() {
            protein.set_is_completely_digested(true);
            match protein.update(self.get_database_connection()) {
                Ok(_) => (),
                Err(err) => summary.log_push(format!("proteomic::models::enzyms::digest_enzym::DigestEnzym::digest(): cannot update protein:\n\t{}", err).as_str())
            }
        }
        return summary;
    }

    fn create_peptide_and_association(&mut self, digest_summary: &mut DigestSummary, protein: &Protein, peptide: &mut Peptide) -> Result<(), String> {
        // create transaction for peptide and association
        let transaction = match self.get_database_connection().transaction() {
            Ok(transaction) => transaction,
            Err(err) => return Err(format!("proteomic::models::enzyms::digest_enzym::DigestEnzym.digest(): Querry error at database_connection.transaction()\n\t{}", handle_postgres_error(&err)))
        };
        transaction.set_rollback(); // set only to commit if no errors occured
        let peptide_created: bool = match peptide.prepared_create(self.get_peptide_create_statement(), self.get_peptide_exists_statement()) {
            Ok(query_ok) => match query_ok {
                QueryOk::Created => true,
                QueryOk::AlreadyExists => false,
                _ => panic!("proteomic::models::enzyms::digest_enzym::DigestEnzym::create_peptide_and_association(): In fact no other QueryOk than QueryOk::Created and QueryOk::AlreadyExists used in Peptide.prepared_create(), this panic should never be reached")
            },
            Err(err) => return Err(format!("proteomic::models::enzyms::digest_enzym::DigestEnzym.create_peptide_and_association(): Querry error at peptide.prepared_create()\n\t{}", err))
        };
        let mut association = PeptideProteinAssociation::new(peptide, protein);
        let association_created: bool = match association.prepared_create(self.get_pp_association_create_statement(), self.get_pp_association_exists_statement()) {
            Ok(query_ok) => match query_ok {
                QueryOk::Created => true,
                QueryOk::AlreadyExists => false,
                _ => panic!("proteomic::models::enzyms::digest_enzym::DigestEnzym::create_peptide_and_association(): In fact no other QueryOk than QueryOk::Created and QueryOk::AlreadyExists used in PeptideProteinAssociation.prepared_create(), this panic should never be reached")
            },
            Err(err) => {
                return Err(format!("proteomic::models::enzyms::digest_enzym::DigestEnzym.create_peptide_and_association(): Querry error at association.prepared_create())\n\t{}", err));
            }
        };
        transaction.set_commit();
        match transaction.finish() {
            Ok(_) => digest_summary.increase_counter(peptide_created, association_created),
            Err(err) => return Err(format!("proteomic::models::enzyms::digest_enzym::DigestEnzym.digest(): Error at transaction commit: {}", handle_postgres_error(&err)))
        }
        return Ok(());
    }

    fn digest_with_hash_set(&self, protein: &mut Protein, aa_sequences: &mut HashSet<String>) {
        /*
         * clone aa_squence and pass it as mutable into replace_all
         * replace every digist_regex-match with with digist_replace (in caseof Trypsin it means add a whitespace between K or T and not P)
         * make the result mutable and split it on whitespaces
         * collect the results as String-vector
         */
        let peptides_without_missed_cleavages: Vec<String> = self.get_digest_regex().split(protein.get_aa_sequence()).map(|peptide| peptide.to_owned()).collect::<Vec<String>>();
        // let mut peptide_position: usize = 0;
        // calculate peptides for missed_cleavages 1 to n + 1 (+1 because explicit boundary)
        'outer: for peptide_idx in 0..peptides_without_missed_cleavages.len() {
            let mut new_peptide_aa_sequence: String = String::new();
            'inner: for number_of_missed_cleavages in 0..(self.get_max_number_of_missed_cleavages() + 1) {
                let temp_idx: usize = peptide_idx + number_of_missed_cleavages as usize;
                if temp_idx < peptides_without_missed_cleavages.len() {
                    new_peptide_aa_sequence.push_str(peptides_without_missed_cleavages.get(temp_idx).unwrap());
                    if self.is_aa_sequence_in_range(&new_peptide_aa_sequence) {
                        aa_sequences.insert(AminoAcid::gerneralize_sequence(&new_peptide_aa_sequence));
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