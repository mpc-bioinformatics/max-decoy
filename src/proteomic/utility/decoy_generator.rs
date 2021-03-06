use std::collections::{HashSet, HashMap};
use std::sync::{Arc, Mutex};
use std::sync::atomic::{AtomicBool, Ordering};
use std::thread;
use std::time;

use time as better_time;

use rand::prelude::*;
use threadpool::ThreadPool;

use proteomic::utility::database_connection::DatabaseConnection;
use proteomic::models::peptides::peptide_interface::PeptideInterface;
use proteomic::models::peptides::peptide::Peptide;
use proteomic::models::peptides::decoy::Decoy;
use proteomic::models::persistable::{Persistable, QueryOk, QueryError};
use proteomic::models::peptides::modified_peptide::{ModifiedPeptide as NewDecoy, PushAminoAcidOk};
use proteomic::models::amino_acids::amino_acid::AminoAcid;
use proteomic::models::amino_acids::amino_acid::AMINO_ACIDS_FOR_DECOY_GENERATION;
use proteomic::models::amino_acids::modification::Modification;


const REPORT_INTERVALL: i64 = 20;   // report all 20 seconds

#[derive(PartialEq)]
pub enum GenerationResult {
    Success,
    Timeout
}

pub struct DecoyGenerator {
    precursor_mass: i64,
    upper_precursor_tolerance_limit: i64,
    lower_precursor_tolerance_limit: i64,
    thread_count: usize,
    max_modifications_per_decoy: u8,
    fixed_modification_map: Arc<HashMap<char, Modification>>,
    variable_modification_map: Arc<HashMap<char, Modification>>,
    one_amino_acid_substitute_map: Arc<HashMap<char, HashMap<char, i64>>>,
    decoys: Arc<Mutex<HashSet<Decoy>>>,
    max_time_for_decoy_generation: i64,
    timeout: Arc<AtomicBool>
}

impl DecoyGenerator {
    pub fn new(precursor_mass: i64, lower_precursor_tolerance_limit: i64, upper_precursor_tolerance_limit: i64, thread_count: usize, max_modifications_per_decoy: u8, fixed_modification_map: &HashMap<char, Modification>, variable_modification_map: &HashMap<char, Modification>, max_time_for_decoy_generation: i64) -> Self {
        return DecoyGenerator{
            precursor_mass: precursor_mass,
            lower_precursor_tolerance_limit: lower_precursor_tolerance_limit,
            upper_precursor_tolerance_limit: upper_precursor_tolerance_limit,
            thread_count: thread_count,
            max_modifications_per_decoy: max_modifications_per_decoy,
            fixed_modification_map: Arc::new(fixed_modification_map.clone()),
            variable_modification_map: Arc::new(variable_modification_map.clone()),
            one_amino_acid_substitute_map: Arc::new(*Self::get_one_amino_acid_substitute_map(fixed_modification_map)),
            decoys: Arc::new(Mutex::new(HashSet::new())),
            max_time_for_decoy_generation: max_time_for_decoy_generation,
            timeout: Arc::new(AtomicBool::new(true))
        }
    }

    pub fn get_lower_precursor_tolerance_limit(&self) -> i64 {
        return self.lower_precursor_tolerance_limit;
    }

    pub fn get_upper_precursor_tolerance_limit(&self) -> i64 {
        return self.upper_precursor_tolerance_limit;
    }

    pub fn get_decoys(&self) -> &Arc<Mutex<HashSet<Decoy>>> {
        return &self.decoys;
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
    fn get_amino_acids_with_weight_less_or_equals_than(weight: i64, fixed_modification_map: &Arc<HashMap<char, Modification>>) -> Box<Vec<char>> {
        let mut amino_acids_and_modification_tupels: Vec<char> = Vec::new();
        for aa_one_letter_code in AMINO_ACIDS_FOR_DECOY_GENERATION {
            let current_amino_acid: AminoAcid = AminoAcid::get(*aa_one_letter_code);
            if let Some(ref modification) = fixed_modification_map.get(aa_one_letter_code) {
                if (current_amino_acid.get_mono_mass() + modification.get_mono_mass()) <= weight {
                    amino_acids_and_modification_tupels.push(*aa_one_letter_code);
                }
            } else {
                if current_amino_acid.get_mono_mass() <= weight {
                    amino_acids_and_modification_tupels.push(*aa_one_letter_code);
                }
            }
        }
        return Box::new(amino_acids_and_modification_tupels);
    }

    pub fn generate_decoys(&self, number_of_decoys_to_generate: usize) -> GenerationResult {
        // set timeout to false before start
        self.timeout.store(false, Ordering::Relaxed);
        // create threadpoll
        let thread_pool = ThreadPool::new(self.thread_count);
        // loop for starting threads
        for thread_id in 0..self.thread_count {
            // create copies of thread safe pointer of DecoyGenerator which can be move into thread
            let fixed_modification_map_ptr = self.fixed_modification_map.clone();
            let variable_modification_map_ptr = self.variable_modification_map.clone();
            let one_amino_acid_substitute_map_ptr = self.one_amino_acid_substitute_map.clone();
            let decoys_ptr = self.decoys.clone();
            let timeout_ptr = self.timeout.clone();
            // copy primitive attributes of DecoyGenerator which can be moved into thread
            let precursor_mass = self.precursor_mass;
            let upper_precursor_tolerance_limit = self.upper_precursor_tolerance_limit;
            let lower_precursor_tolerance_limit = self.lower_precursor_tolerance_limit;
            let max_modifications_per_decoy = self.max_modifications_per_decoy;
            // start thread
            thread_pool.execute(move||{
                let conn: postgres::Connection = DatabaseConnection::get_database_connection();
                // create random number generator
                let mut rng = rand::thread_rng();
                // endless loop with label 'decoy_loop
                'decoy_loop: loop {
                    match decoys_ptr.lock() {
                        Ok(decoys) => if decoys.len() >= number_of_decoys_to_generate { break 'decoy_loop; },
                        Err(_) => panic!("proteomic::utility::decoy_generator::DecoyGenerator.generate_decoys(): try to lock poisened mutex for decoys")
                    };
                    if timeout_ptr.load(Ordering::Relaxed) { break 'decoy_loop; }
                    // create new empty decoy
                    let mut new_decoy: NewDecoy = NewDecoy::new_decoy(precursor_mass, lower_precursor_tolerance_limit, upper_precursor_tolerance_limit);
                    // let distribution_array = *Self::generate_amino_acid_distribution_array();
                    // repeat until new_decoy's weight greate then upper weight limit
                    'amino_acid_loop: loop {
                        // pick amino acids one letter code at index
                        //let aa_one_letter_code: char = *distribution_array.choose(&mut rng).unwrap();
                        let aa_one_letter_code: char = *AMINO_ACIDS_FOR_DECOY_GENERATION.choose(&mut rng).unwrap();
                        // one letter code to amino acid
                        let random_amino_acid: AminoAcid = AminoAcid::get(aa_one_letter_code);
                        let modification_option = match fixed_modification_map_ptr.get(&random_amino_acid.get_one_letter_code()){
                            Some(ref modification) => Some((*modification).clone()),
                            None => None
                        };
                        match new_decoy.push_amino_acid_and_fix_modification(&random_amino_acid, &modification_option) {
                            Ok(push_ok) => match push_ok {
                                PushAminoAcidOk::GreaterThenMassTolerance => break 'amino_acid_loop,
                                _ => ()
                            },
                            Err(push_err) => panic!("proteomic::utility::decoy_generator::DecoyGenerator.generate_decoys(): Error at new_decoy.push_amino_acid_and_fix_modification: {}", push_err)
                        }
                    }
                    if new_decoy.hits_mass_tolerance() {
                        let mut decoy = new_decoy.to_decoy();
                        decoy.set_modification_summary(new_decoy.get_modification_summary_for_header().as_str());
                        if Self::save_new_decoy(&conn, &thread_id, &mut decoy) {
                            match decoys_ptr.lock() {
                                Ok(mut decoys) => {
                                    decoys.insert(decoy);
                                    continue 'decoy_loop;
                                },
                                Err(_) => panic!("proteomic::utility::decoy_generator::DecoyGenerator.generate_decoys(): try to lock poisoned mutex for decoys")
                            }
                        }
                    }
                    if new_decoy.swap_amino_acids_to_hit_mass_tolerance(one_amino_acid_substitute_map_ptr.as_ref(), fixed_modification_map_ptr.as_ref(), max_modifications_per_decoy, variable_modification_map_ptr.as_ref())
                    {
                        let mut decoy = new_decoy.to_decoy();
                        decoy.set_modification_summary(new_decoy.get_modification_summary_for_header().as_str());
                        if Self::save_new_decoy(&conn, &thread_id, &mut decoy) {
                            match decoys_ptr.lock() {
                                Ok(mut decoys) => {
                                    decoys.insert(decoy);
                                    continue 'decoy_loop;
                                },
                                Err(_) => panic!("proteomic::utility::decoy_generator::DecoyGenerator.generate_decoys(): try to lock poisoned mutex for decoys")
                            }
                        }
                    }
                }
            });
        }
        let start_at_sec = better_time::now().to_timespec().sec;                // lock start sec
        let stop_at_sec = start_at_sec + self.max_time_for_decoy_generation;    // lock stop sec
        let mut next_report_at = start_at_sec + REPORT_INTERVALL;               // next second to report
        thread::sleep(time::Duration::from_secs(2));                            // wait for threads to start
        while thread_pool.active_count() > 0 {
            let now_in_seconds = better_time::now().to_timespec().sec;          // current second
            if now_in_seconds >= next_report_at {
                next_report_at = now_in_seconds + REPORT_INTERVALL;
                match self.decoys.lock() {
                    Ok(decoys) => {
                        println!("decoy generation progress: {} decoys generated", decoys.len());
                    },
                    Err(_) => panic!("proteomic::utility::decoy_generator::DecoyGenerator.generate_decoys(): Tried to lock poisoned mutex for decoys")
                }
            }
            // signal threads to stop and break the loop if max generation time is reached
            if now_in_seconds >= stop_at_sec {
                self.timeout.store(true, Ordering::Relaxed);
                println!("proteomic::utility::decoy_generator::DecoyGenerator.generate_decoys(): Running out of time, stop threads");
                break;
            }
        }
        // wait for generator threads to stop
        thread_pool.join();
        if !self.timeout.load(Ordering::Relaxed) {
            return GenerationResult::Success;
        } else {
            return GenerationResult::Timeout
        }
    }

    fn save_new_decoy(conn: &postgres::Connection, thread_id: &usize, decoy: &mut Decoy) -> bool {
        let mut error_occured = false;
        if !decoy.is_peptide(&conn) {
            // try 3 times to create the decoy
            for _ in 0..3 {
                match decoy.create(&conn) {
                    Ok(query_ok) => {
                        let result = match query_ok {
                            // return true only if decoy is newly created. decoys which not newly created are already added from other threads or
                            // add during previous target/decoy search
                            QueryOk::Created => true,
                            _ => false
                        };
                        // if previous
                        if error_occured {
                            println!("proteomic::utility::decoy_generator::DecoyGenerator.generate_decoys()_thread_{}: Previous error resolved", thread_id);
                        }
                        return result;
                    }
                    Err(query_err) => match query_err {
                        // NoReturn is mostly a indication for a conflict between two or more threads which try to create the same decoy
                        QueryError::NoReturn => {
                            error_occured = true;
                            // wait a random time between 1 and 6 seconds so hopefully one thread get the chance to create the decoy before the other threads
                            let mut rng = rand::thread_rng();
                            let wait_for: u64 = rng.gen_range(1, 6);
                            println!("proteomic::utility::decoy_generator::DecoyGenerator.generate_decoys()_thread_{}: Conflict with other threads which try to create the same decoy. Try again in {} seconds", thread_id, wait_for);
                            thread::sleep(time::Duration::from_secs(wait_for));
                            continue;
                        },
                        _ => {
                            println!("proteomic::utility::decoy_generator::DecoyGenerator.generate_decoys()_thread_{}: Unresoveable error occured: {}, continue with next decoy", thread_id, query_err);
                            return false;
                        }
                    }
                }
            }
        }
        if error_occured {
            println!("proteomic::utility::decoy_generator::DecoyGenerator.generate_decoys()_thread_{}: Previous error could not resolved, continue with next decoy", thread_id);
        }
        return false;
    }

    pub fn vary_targets(&self, targets: &Vec<Peptide>, number_of_decoys_to_generate: usize) -> Box<HashSet<Decoy>> {
        let mut decoys: Box<HashSet<Decoy>> = Box::new(HashSet::new());
        let conn: postgres::Connection = DatabaseConnection::get_database_connection();
        let mut rng = rand::thread_rng();
        'targets_loop: for target in targets.iter() {
            let mut aa_sequence: Vec<char> = target.get_aa_sequence().chars().collect();
            'shuffle_loop: for _ in 0..1000 {
                aa_sequence.shuffle(&mut rng);
                let aa_sequence_as_string: String = aa_sequence.iter().map(|amino_acid| *amino_acid).collect::<String>();
                match Peptide::exists_where(&conn, "aa_sequence = $1", &[&aa_sequence_as_string]) {
                    Ok(_) => continue 'shuffle_loop,  // if decoy is peptide, continue with next iteration
                    Err(err) => match err {
                        QueryError::NoMatch => (),  // if no match is found continue with this decoy
                        _ => panic!("proteomic::utility::decoy_generator::DecoyGenerator.vary_decoy() could not check if shuffled target is peptide: {}", err)
                    }
                }
                let mut new_decoy: NewDecoy = NewDecoy::decoy_from_string(aa_sequence_as_string.as_str(), self.precursor_mass, self.get_lower_precursor_tolerance_limit(), self.get_upper_precursor_tolerance_limit(), self.fixed_modification_map.as_ref());
                if new_decoy.hits_mass_tolerance() {
                    let mut decoy = new_decoy.to_decoy();
                    if !decoy.is_peptide(&conn) {
                        match decoy.create(&conn) {
                            Ok(_) => { decoys.insert(decoy); },    // wrap insert into block, to 'suppress' return
                            Err(err) => println!("proteomic::utility::decoy_generator::DecoyGenerator.vary_targets(): Could not create decoy: {}", err)
                        }
                    }
                }
                // at this point, we could also swap amino acids and try to apply modifications
                if decoys.len() >= number_of_decoys_to_generate { break 'shuffle_loop; };
            }
        }
        return decoys;
    }

    /// calculates a substitution map for swapping amino acids with each other.
    /// so one can lookup which difference in weight a substitution of amino acid x with y has.
    /// this function also considers fixed modifications
    pub fn get_one_amino_acid_substitute_map(fixed_modification_map: &HashMap<char, Modification>) -> Box<HashMap<char, HashMap<char, i64>>> {
        let mut substitution_map: HashMap<char, HashMap<char, i64>> = HashMap::new();
        for aa_origin in AMINO_ACIDS_FOR_DECOY_GENERATION.iter() {
            let mut differences_in_weight: HashMap<char, i64> = HashMap::new();
            let mut aa_origin_weight: i64 = AminoAcid::get(*aa_origin).get_mono_mass();
            aa_origin_weight += match fixed_modification_map.get(aa_origin){
                Some(ref modification) => modification.get_mono_mass(),
                None => 0
            };
            for aa_replacement in AMINO_ACIDS_FOR_DECOY_GENERATION.iter() {
                if *aa_replacement == *aa_origin {
                    continue;
                }
                let mut aa_replacement_weight: i64 = AminoAcid::get(*aa_replacement).get_mono_mass();
                aa_replacement_weight += match fixed_modification_map.get(aa_replacement){
                    Some(ref modification) => modification.get_mono_mass(),
                    None => 0
                };
                differences_in_weight.insert(*aa_replacement, aa_replacement_weight - aa_origin_weight);
            }
            substitution_map.insert(*aa_origin, differences_in_weight);
        }
        return Box::new(substitution_map);
    }
}