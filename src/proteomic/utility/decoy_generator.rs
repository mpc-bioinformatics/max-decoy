extern crate rand;
extern crate threadpool;


use std::collections::HashMap;
use std::sync::Mutex;
use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering};

use self::rand::Rng;
use self::rand::seq::SliceRandom;
use self::threadpool::ThreadPool;

use proteomic::models::mass;
// use proteomic::utility::mass::NeutralLoss;
// use proteomic::models::peptide::Peptide;
use proteomic::models::decoys::decoy::Decoy;
use proteomic::models::decoys::new_decoy::{NewDecoy, NewDecoyError};
use proteomic::models::decoys::base_decoy::BaseDecoy;
use proteomic::models::amino_acids::amino_acid::AminoAcid;
use proteomic::models::amino_acids::modification::Modification;

const AMINO_ACIDS_FOR_DECOY_GENERATION: &'static [char] = &['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'J', 'K', 'M', 'F', 'P', 'O', 'S', 'T', 'V', 'W', 'Y'];

pub struct DecoyGenerator {
    upper_weight_limit: i64,
    lower_weight_limit: i64,
    decoy_limit: usize,
    thread_count: usize,
    decoy_counter: Arc<AtomicUsize>,
    fixed_modification_map: Arc<Mutex<HashMap<char, Modification>>>,
    variable_modification_map: Arc<Mutex<HashMap<char, Modification>>>,
}

impl DecoyGenerator {
    pub fn new(weight: f64, upper_limit_ppm: i64, lower_limit_ppm: i64, decoy_limit: usize, thread_count: usize, fixed_modification_map: &HashMap<char, Modification>, variable_modification_map: &HashMap<char, Modification>) -> Self {
        let weight_as_int = mass::convert_mass_to_int(weight);
        let weight_upper_limit_ppm = weight_as_int / 1000000 * upper_limit_ppm as i64;
        let weight_lower_limit_ppm = weight_as_int / 1000000 * lower_limit_ppm as i64;
        return DecoyGenerator{
            upper_weight_limit: weight_as_int + weight_upper_limit_ppm,
            lower_weight_limit: weight_as_int - weight_lower_limit_ppm,
            decoy_limit: decoy_limit,
            thread_count: thread_count,
            decoy_counter: Arc::new(AtomicUsize::new(0)),
            fixed_modification_map: Arc::new(Mutex::new(fixed_modification_map.clone())),
            variable_modification_map: Arc::new(Mutex::new(variable_modification_map.clone()))
        }
    }

    // generate array which holds natural distribution of amino acids (published by UniProt)
    fn generate_amino_acid_distribution_array() -> Box<Vec<char>> {
        // create empty array
        let mut distribution_array: Vec<char> = Vec::new();
        // loop through all amino acids which are allowed in decoy generation
        for aa_one_letter_code in AMINO_ACIDS_FOR_DECOY_GENERATION {
            let current_amino_acid: AminoAcid = AminoAcid::get(*aa_one_letter_code);
            for _ in 0..current_amino_acid.get_distribution() {
                distribution_array.push(*aa_one_letter_code);
            }
        }
        // return boxed array
        return Box::new(distribution_array);
    }

    // generate array with amino acids which have a mass less then a specific weight
    // by considering the fixed modifications
    fn get_amino_acids_which_weight_less_than(weight: i64, fixed_modification_ptr: &Arc<Mutex<HashMap<char, Modification>>>) -> Box<Vec<char>> {
        let mut amino_acids: Vec<char> = Vec::new();
        for aa_one_letter_code in AMINO_ACIDS_FOR_DECOY_GENERATION {
            let current_amino_acid: AminoAcid = AminoAcid::get(*aa_one_letter_code);
            match fixed_modification_ptr.lock() {
                Ok(modifications) => {
                    match modifications.get(aa_one_letter_code) {
                        Some(modification) => {
                            if (current_amino_acid.get_mono_mass() + modification.get_mono_mass()) <= weight {
                                amino_acids.push(current_amino_acid.get_one_letter_code());
                            }
                        },
                        None => {
                            if current_amino_acid.get_mono_mass() <= weight {
                                amino_acids.push(current_amino_acid.get_one_letter_code());
                            }
                        }
                    }
                }
                Err(_) => println!("ERROR DecoyGenerator: Could not gather lock for fixed modifications")
            }
        }
        return Box::new(amino_acids);
    }

    pub fn generate_decoys(&self) {
        // create threadpoll
        let thread_pool = ThreadPool::new(self.thread_count);
        // loop for starting threads
        for _ in 0..self.thread_count {
            // create copies of thread safe pointer of DecoyGenerator which can be move into thread
            let decoy_counter_ptr = self.decoy_counter.clone();
            let fixed_modification_map_ptr = self.fixed_modification_map.clone();
            let variable_modification_map_ptr = self.variable_modification_map.clone();
            // copy primitive attributes of DecoyGenerator which can be moved into thread
            let upper_weight_limit = self.upper_weight_limit;
            let lower_weight_limit = self.lower_weight_limit;
            let decoy_limit = self.decoy_limit;
            // start thread
            thread_pool.execute(move||{
                // create random number generator
                let mut rng = rand::thread_rng();
                // get haviest and lightest amino acids
                let haviest_amino_acid: AminoAcid = AminoAcid::get_haviest();
                let mut haviest_amino_acid_mass_with_modification: i64 = haviest_amino_acid.get_mono_mass();
                match fixed_modification_map_ptr.lock() {
                    Ok(modifications) => {
                        match modifications.get(&haviest_amino_acid.get_one_letter_code()){
                            Some(modification) => {
                                haviest_amino_acid_mass_with_modification += modification.get_mono_mass();
                            },
                            None => {
                                haviest_amino_acid_mass_with_modification += 0;
                            }
                        }
                    },
                    Err(_) => println!("ERROR DecoyGenerator: Could not gather lock for fixed modifications")
                }
                // let lightest_amino_acid = AminoAcid::get_lightest();
                // endless loop with label 'decoy_loop
                'decoy_loop: loop {
                    // create new empty decoy
                    let mut new_decoy: NewDecoy = NewDecoy::new(lower_weight_limit, upper_weight_limit);
                    // get array of fitting amino acids
                    let distribution_array = *Self::generate_amino_acid_distribution_array();
                    // repeat until new_decoy's weight is smaller than the haviest amino acid. so
                    while new_decoy.get_weight() >= haviest_amino_acid_mass_with_modification {
                        // pick amino acids one letter code at index
                        let aa_one_letter_code: char = *distribution_array.choose(&mut rng).unwrap();
                        // one letter code to amino acid
                        let random_amino_acid: AminoAcid = AminoAcid::get(aa_one_letter_code);
                        // gather lock on fixed modification
                        match fixed_modification_map_ptr.lock() {
                            Ok(fixed_modification_map) => {
                                // add amino acid and modification to decoy
                                match new_decoy.push_amino_acid_with_modification(
                                    &random_amino_acid,
                                    &fixed_modification_map.get(&random_amino_acid.get_one_letter_code())
                                ) {
                                    Ok(_) => (),
                                    Err(err) => {
                                        match err {
                                            NewDecoyError::ModificationDoesNotMatchToAminoAcid => panic!("tried to push '{}' with an incompatible modification", random_amino_acid.get_one_letter_code())
                                        }
                                    }
                                }
                                return ();
                            },
                            // at this point something is really faulty, .lock() should wait until lock is free
                            // so if you get here the lock ist maybe poisoned
                            Err(_) => println!("ERROR DecoyGenerator: Could not gather lock for fixed modifications")
                        }
                        println!("{} => {}", new_decoy.get_aa_sequence(), new_decoy.get_weight());
                    }
                    // fill up with amino acids which fits the remaining weights
                    let mut fitting_amino_acids: Vec<char> = *Self::get_amino_acids_which_weight_less_than(new_decoy.get_weight(), &fixed_modification_map_ptr);
                    while fitting_amino_acids.len() > 0 {
                        let aa_one_letter_code: char = *distribution_array.choose(&mut rng).unwrap();
                        let random_amino_acid: AminoAcid = AminoAcid::get(aa_one_letter_code);
                        match fixed_modification_map_ptr.lock() {
                            Ok(fixed_modification_map) => {
                                match new_decoy.push_amino_acid_with_modification(
                                    &random_amino_acid,
                                    &fixed_modification_map.get(&random_amino_acid.get_one_letter_code())
                                ) {
                                    Ok(_) => (),
                                    Err(err) => {
                                        match err {
                                            NewDecoyError::ModificationDoesNotMatchToAminoAcid => panic!("tried to push '{}' with an incompatible modification", random_amino_acid.get_one_letter_code())
                                        }
                                    }
                                }
                                return ();
                            },
                            // at this point something is really faulty, .lock() should wait until lock is free
                            // so if you get here the lock ist maybe poisoned
                            Err(_) => println!("ERROR DecoyGenerator: Could not gather lock for fixed modifications")
                        }
                        fitting_amino_acids = *Self::get_amino_acids_which_weight_less_than(new_decoy.get_weight(), &fixed_modification_map_ptr);
                    }
                    if new_decoy.hits_mass_tolerance() {
                        decoy_counter_ptr.fetch_add(1, Ordering::Relaxed);
                        let base_decoy: BaseDecoy = new_decoy.as_base_decoy();
                        println!("{}", base_decoy.to_string());
                        println!("{}", new_decoy.as_modified_decoy(&base_decoy).to_string());
                        if decoy_counter_ptr.load(Ordering::Relaxed) >= 1000 {
                            break 'decoy_loop;
                        }
                    } else {
                        println!("drop decoy with weight {} ({})", new_decoy.get_weight(), new_decoy.get_modified_weight());
                    }
                }
            });
        }
        println!("tp.count => {}", thread_pool.active_count());
        thread_pool.join();
    }

    fn is_decoy_in_range(&self, decoy: &Decoy) -> bool {
        return (self.lower_weight_limit <= decoy.get_weight()) & (decoy.get_weight() <= self.upper_weight_limit)
    }

    fn is_weight_in_range(&self, weight: i64) -> bool {
        return (self.lower_weight_limit <= weight) & (weight <= self.upper_weight_limit)
    }

    fn find_amino_acid_which_fits(&self, weight: i64) -> Option<char> {
        for one_letter_code in AMINO_ACIDS_FOR_DECOY_GENERATION {
            let current_amino_acid = AminoAcid::get(*one_letter_code);
            if self.is_weight_in_range(weight + current_amino_acid.get_mono_mass()) {
                return Some(*one_letter_code);
            }
        }
        return None;
    }
}