extern crate regex;

use std::collections::HashSet;
use std::iter::FromIterator;

pub struct Trypsin {
    digist_regex: regex::Regex,
    digest_replace: &'static str,
    missed_cleavages: usize,
    min_peptide_length: usize,
    max_peptide_length: usize
}

pub trait DigstEnzym {
    fn new(missed_cleavages: usize, min_peptide_length: usize, max_peptide_length: usize) -> Self;

    fn digest(&self, protein: &mut super::protein::Protein) -> HashSet<String>;


    
    fn digest_internal(&self, peptides: &Vec<String>, split_positions: &Vec<usize>, number_of_missed_cleavages: usize, incl_arr: &mut Vec<bool>, index: usize, mut current_peptide_set: HashSet<String>) -> HashSet<String> {
        if split_positions.len() < number_of_missed_cleavages + index { return current_peptide_set; }
        if number_of_missed_cleavages == 0 {
            let it = split_positions.iter().zip(incl_arr.iter()).filter_map(|(val, incl)|
                if *incl { Some(val) } else { None }
            );
            let mut missed_cleavage_combination: Vec<usize> = Vec::new();
            for val in it { 
                missed_cleavage_combination.push(*val);
            }
            return self.process_combinations(peptides, &missed_cleavage_combination, current_peptide_set);
        }
        
        incl_arr[index] = true;
        current_peptide_set = self.digest_internal(peptides, split_positions, number_of_missed_cleavages-1, incl_arr, index+1, current_peptide_set);
        incl_arr[index] = false;
        
        current_peptide_set = self.digest_internal(peptides, split_positions, number_of_missed_cleavages, incl_arr, index+1, current_peptide_set);

        return current_peptide_set;
    }

    fn process_combinations(&self, peptides: &Vec<String>, missed_cleavage_combination: &Vec<usize>, mut current_peptide_set: HashSet<String>) -> HashSet<String>{
        let mut new_peptide = String::new();
        for i in 0..missed_cleavage_combination.len() {
            new_peptide.push_str(peptides[missed_cleavage_combination[i]].clone().as_mut_str());
            if (i + 1) < (missed_cleavage_combination.len()) && missed_cleavage_combination[i+1] != missed_cleavage_combination[i] + 1 {
                current_peptide_set.insert(new_peptide);
                new_peptide = String::new();
            }
        }
        current_peptide_set.insert(new_peptide);
        return current_peptide_set;
    }
}


impl DigstEnzym for Trypsin {
    fn new (missed_cleavages: usize, min_peptide_length: usize, max_peptide_length: usize) -> Trypsin {
        Trypsin {
            digist_regex: regex::Regex::new(r"(?P<before>(K|R))(?P<after>[^P])").unwrap(),
            digest_replace: "$before $after",
            missed_cleavages: missed_cleavages,
            min_peptide_length: min_peptide_length,
            max_peptide_length: max_peptide_length
        }
    }

    fn digest(&self, protein: &mut super::protein::Protein) -> HashSet<String> {
        /*
         * clone aa_squence and pass it as mutable into replace_all
         * replace every digist_regex-match with with digist_replace (in caseof Trypsin it means add a whitespace between K or T and not P)
         * make the result mutable and split it on whitespaces
         * collect the results as String-vector
         */
        let peptides_without_missed_cleavages: Vec<String> = self.digist_regex.replace_all(protein.get_aa_sequence().clone().as_mut_str(), self.digest_replace).to_mut().split(" ").map(|peptide| peptide.to_owned()).collect::<Vec<String>>();
        let mut peptides: HashSet<String> = HashSet::from_iter(peptides_without_missed_cleavages.clone());
        // reccuring variables for calculation of missed-cleavages-combinations
        let mut incl_arr: Vec<bool> = vec![false; peptides.len()];
        let split_positions: Vec<usize> = (0..peptides.len()).collect();
        // calculate peptides for missed_cleavages 1 to n + 1 (+1 because explicit boundary)
        for number_of_missed_cleavages in 1..(self.missed_cleavages + 1) {
            peptides = self.digest_internal(&peptides_without_missed_cleavages, &split_positions, number_of_missed_cleavages, &mut incl_arr, 0, peptides);
        }

        for pep in peptides {
            println!("{}", pep);
        }
        return HashSet::new();
    }
}