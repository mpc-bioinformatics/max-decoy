extern crate regex;

use std::collections::HashSet;

use proteomic::models::protein::Protein;
use proteomic::models::peptide::Peptide;

pub struct Trypsin {
    name: String,
    digist_regex: regex::Regex,
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
    fn get_digist_regex(&self) -> &regex::Regex;
    fn get_digist_replace(&self) -> &'static str;
    fn get_min_peptide_length(&self) -> usize;
    fn get_max_peptide_length(&self) -> usize;

    fn digest(&self, protein: &mut Protein) -> HashSet<Peptide> {
        /*
         * clone aa_squence and pass it as mutable into replace_all
         * replace every digist_regex-match with with digist_replace (in caseof Trypsin it means add a whitespace between K or T and not P)
         * make the result mutable and split it on whitespaces
         * collect the results as String-vector
         */
        let peptides_without_missed_cleavages: Vec<String> = self.get_digist_regex().replace_all(protein.get_aa_sequence().clone().as_mut_str(), self.get_digist_replace()).to_mut().split(" ").map(|peptide| peptide.to_owned()).collect::<Vec<String>>();
        let mut peptides: HashSet<Peptide> = HashSet::new();
        for peptide_aa_sequence in &peptides_without_missed_cleavages {
            if self.is_peptide_aa_sequence_in_range(&peptide_aa_sequence){
                peptides.insert(Peptide::new(peptide_aa_sequence.to_owned(), 0, String::from(self.get_name())));
            }
        }
        // reccuring variables for calculation of missed-cleavages-combinations
        let mut incl_arr: Vec<bool> = vec![false; peptides.len()];
        let cleavage_positions: Vec<usize> = (0..peptides.len()).collect();
        // calculate peptides for missed_cleavages 1 to n + 1 (+1 because explicit boundary)
        for number_of_missed_cleavages in 1..(self.get_max_number_of_missed_cleavages() + 1) {
            self.digest_internal(&peptides_without_missed_cleavages, &cleavage_positions, number_of_missed_cleavages, &mut incl_arr, 0, &mut peptides);
        }
        return peptides;
    }

    /*
     * calculates combination n choose k
     * n = cleavage_positions
     * k = number_of_missed_cleavages
     * source: https://rosettacode.org/wiki/Combinations#Rust
     * modification: does not print combination it passes the combination to process combination instead
     */
    fn digest_internal(&self, peptides: &Vec<String>, cleavage_positions: &Vec<usize>, number_of_missed_cleavages: usize, incl_arr: &mut Vec<bool>, index: usize, current_peptide_set: &mut HashSet<Peptide>) {
        if cleavage_positions.len() < number_of_missed_cleavages + index { return; }
        if number_of_missed_cleavages == 0 {
            let it = cleavage_positions.iter().zip(incl_arr.iter()).filter_map(|(val, incl)|
                if *incl { Some(val) } else { None }
            );
            let mut missed_cleavage_combination: Vec<usize> = Vec::new();
            for val in it { 
                missed_cleavage_combination.push(*val);
            }
            self.process_combinations(peptides, &missed_cleavage_combination, current_peptide_set);
            return;
        }
        
        incl_arr[index] = true;
        self.digest_internal(peptides, cleavage_positions, number_of_missed_cleavages-1, incl_arr, index+1, current_peptide_set);
        incl_arr[index] = false;
        
        self.digest_internal(peptides, cleavage_positions, number_of_missed_cleavages, incl_arr, index+1, current_peptide_set);

        return;
    }

    /*
     * adds the peptides in missed_cleavage_combination to current_peptide_set. the values in missed_cleavages_combination represents the missed cleavage positions
     * lets assume missed cleavages is something like [1,4,7], the function will add peptide on position 1, 4 and 7 to current_peptide_set.
     * if two or more indices in missed_cleavage_combination are sequent, the function will glue the peptides together.
     * lets assume missed cleavages is something like [1,2,6], the function will glue the peptides 1, 2 and 3 together and add the contatenation to current_peptide_set also peptide 6 concatenated with 7 will added.
     */
    fn process_combinations(&self, peptides: &Vec<String>, missed_cleavage_combination: &Vec<usize>, current_peptide_set: &mut HashSet<Peptide>){
        let mut new_peptide_aa_sequence = String::new();
        for i in 0..missed_cleavage_combination.len() {
            // if new_peptide_aa_sequence has length 0, add both peptides around the cleavage position to the new_peptide_aa_sequence (glue them together)
            if new_peptide_aa_sequence.len() == 0 {
                new_peptide_aa_sequence.push_str(peptides[missed_cleavage_combination[i]].clone().as_mut_str());
                new_peptide_aa_sequence.push_str(peptides[missed_cleavage_combination[i] + 1].clone().as_mut_str());
            } 
            // add only the peptide on the left of the cleavage position to the peptide, because the right one is added in the iteration step before
            else {
                new_peptide_aa_sequence.push_str(peptides[missed_cleavage_combination[i] + 1].clone().as_mut_str());
            }
            /*
             * if i + 1 is less then the length of missed_cleavage_combination and the next missed cleavages position is NOT sequent to the current
             * add the new_peptide_aa_sequence to the results and begin with a new one
             */
            if (i + 1) < (missed_cleavage_combination.len()) && missed_cleavage_combination[i+1] != missed_cleavage_combination[i] + 1 {
                // println!("{} => {}", new_peptide_aa_sequence, self.is_peptide_in_range(&new_peptide_aa_sequence));
                if self.is_peptide_aa_sequence_in_range(&new_peptide_aa_sequence) {
                    current_peptide_set.insert(Peptide::new(new_peptide_aa_sequence, 0, String::from(self.get_name())));
                }
                new_peptide_aa_sequence = String::new();
            }
        }
        // println!("{} => {}", new_peptide_aa_sequence, self.is_peptide_in_range(&new_peptide_aa_sequence));
        if self.is_peptide_aa_sequence_in_range(&new_peptide_aa_sequence) {
            current_peptide_set.insert(Peptide::new(new_peptide_aa_sequence, 0, String::from(self.get_name())));
        }
    }

    /*
     * checks if the passed peptide has the correct length
     * note: maybe this method ca be beatufied in the feature with https://doc.rust-lang.org/std/ops/struct.Range.html#method.contains
     */
    fn is_peptide_aa_sequence_in_range(&self, new_peptide_aa_sequence: &String) -> bool {
        return self.get_min_peptide_length() <= new_peptide_aa_sequence.len() && new_peptide_aa_sequence.len() <= self.get_max_peptide_length()
    }
}


impl DigestEnzym for Trypsin {
    fn new (max_number_of_missed_cleavages: usize, min_peptide_length: usize, max_peptide_length: usize) -> Trypsin {
        Trypsin {
            name: String::from("Trypsin"),
            digist_regex: regex::Regex::new(r"(?P<before>(K|R))(?P<after>[^P])").unwrap(),
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

    fn get_digist_regex(&self) -> &regex::Regex {
        return &self.digist_regex;
    }

    fn get_digist_replace(&self) -> &'static str {
        return self.digest_replace;
    }

    fn get_min_peptide_length(&self) -> usize {
        return self.min_peptide_length;
    }

    fn get_max_peptide_length(&self) -> usize {
        return self.max_peptide_length;
    }
}