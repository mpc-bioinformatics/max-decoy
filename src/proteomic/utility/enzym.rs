extern crate onig;
extern crate postgres;

use std::collections::HashSet;

use proteomic::models::protein::Protein;
use proteomic::models::peptide::Peptide;
use proteomic::models::peptide_protein_association::PeptideProteinAssociation;

pub struct Trypsin {
    name: String,
    digist_regex: onig::Regex,
    digest_replace: &'static str,
    max_number_of_missed_cleavages: usize,
    // replace this with a range in the feature: https://doc.rust-lang.org/std/ops/struct.Range.html#method.contains
    min_peptide_length: usize,
    max_peptide_length: usize
}

pub trait DigestEnzym {
    fn new(max_number_of_missed_cleavages: usize, min_peptide_length: usize, max_peptide_length: usize) -> Self;
    fn get_name(&self) -> &str;
    fn get_max_number_of_missed_cleavages(&self) -> usize;
    fn get_digest_regex(&self) -> &onig::Regex;
    fn get_digest_replace(&self) -> &'static str;
    fn get_min_peptide_length(&self) -> usize;
    fn get_max_peptide_length(&self) -> usize;

    fn digest(&self, protein: &mut Protein, peptides: &mut HashSet<Peptide>, database_connection: &postgres::Connection) {
        /*
         * clone aa_squence and pass it as mutable into replace_all
         * replace every digist_regex-match with with digist_replace (in caseof Trypsin it means add a whitespace between K or T and not P)
         * make the result mutable and split it on whitespaces
         * collect the results as String-vector
         */
        let peptides_without_missed_cleavages: Vec<String> = self.get_digest_regex().split(protein.get_aa_sequence().clone().as_mut_str()).map(|peptide| peptide.to_owned()).collect::<Vec<String>>();
        let mut peptide_position: usize = 0;
        // calculate peptides for missed_cleavages 1 to n + 1 (+1 because explicit boundary)
        'outer: for peptide_idx in 0..peptides_without_missed_cleavages.len() {
            let mut new_peptide_aa_sequence: String = String::new();
            'inner: for number_of_missed_cleavages in 0..(self.get_max_number_of_missed_cleavages() + 1) {
                let temp_idx: usize = peptide_idx + number_of_missed_cleavages;
                if temp_idx < peptides_without_missed_cleavages.len() {
                    new_peptide_aa_sequence.push_str(peptides_without_missed_cleavages.get(temp_idx).unwrap());
                    if self.is_aa_sequence_in_range(&new_peptide_aa_sequence) {
                        let mut peptide = Peptide::new(new_peptide_aa_sequence.clone(), self.get_name().to_owned(), number_of_missed_cleavages as i32);
                        peptide.save(&database_connection);
                        peptides.insert(peptide);
                    }
                } else {
                    break;
                }
            }
            peptide_position += peptides_without_missed_cleavages[peptide_idx].len();
        }
    }

    fn is_aa_sequence_in_range(&self, aa_sequence: &String) -> bool {
        return self.get_min_peptide_length() <= aa_sequence.len() && aa_sequence.len() <= self.get_max_peptide_length();
    }
}


impl DigestEnzym for Trypsin {
    fn new (max_number_of_missed_cleavages: usize, min_peptide_length: usize, max_peptide_length: usize) -> Trypsin {
        Trypsin {
            name: String::from("Trypsin"),
            digist_regex: onig::Regex::new(r"(?<=[KR])(?!P)").unwrap(),
            digest_replace: "$before $after",
            max_number_of_missed_cleavages: max_number_of_missed_cleavages,
            min_peptide_length: min_peptide_length,
            max_peptide_length: max_peptide_length
        }
    }

    fn get_name(&self) -> &str {
        return self.name.as_str();
    }

    fn get_max_number_of_missed_cleavages(&self) -> usize{
        return self.max_number_of_missed_cleavages;
    }

    fn get_digest_regex(&self) -> &onig::Regex {
        return &self.digist_regex;
    }

    fn get_digest_replace(&self) -> &'static str {
        return self.digest_replace;
    }

    fn get_min_peptide_length(&self) -> usize {
        return self.min_peptide_length;
    }

    fn get_max_peptide_length(&self) -> usize {
        return self.max_peptide_length;
    }
}