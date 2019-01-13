extern crate rand;
extern crate threadpool;


use std::collections::{HashSet, HashMap};
use std::sync::Mutex;
use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::thread;

use self::rand::seq::SliceRandom;
use self::threadpool::ThreadPool;

use proteomic::models::mass;
use proteomic::utility::database_connection::DatabaseConnection;
use proteomic::utility::combinations::n_choose_k::NChooseK;
// use proteomic::utility::mass::NeutralLoss;
use proteomic::models::peptide::Peptide;
use proteomic::models::persistable::{Persistable, QueryOk, QueryError};
use proteomic::models::decoys::decoy::{Decoy, PlainDecoy};
use proteomic::models::decoys::new_decoy::{NewDecoy, NewDecoyError};
use proteomic::models::decoys::base_decoy::BaseDecoy;
use proteomic::models::amino_acids::amino_acid::AminoAcid;
use proteomic::models::amino_acids::modification::Modification;

const AMINO_ACIDS_FOR_DECOY_GENERATION: &'static [char] = &['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'J', 'K', 'M', 'F', 'P', 'O', 'S', 'T', 'V', 'W', 'Y'];

pub struct DecoyGenerator {
    upper_weight_limit: i64,
    lower_weight_limit: i64,
    thread_count: usize,
    decoy_counter: Arc<AtomicUsize>,
    max_modifications_per_decoy: i32,
    fixed_modification_map: Arc<HashMap<char, Modification>>,
    variable_modification_map: Arc<HashMap<char, Modification>>,
    one_amino_acid_substitute_map: Arc<HashMap<char, HashMap<char, i64>>>
}

impl DecoyGenerator {
    pub fn new(weight: f64, lower_mass_limit_ppm: i64, upper_mass_limit_ppm: i64, thread_count: usize, max_modifications_per_decoy: i32, fixed_modification_map: &HashMap<char, Modification>, variable_modification_map: &HashMap<char, Modification>) -> Self {
        return DecoyGenerator{
            upper_weight_limit: mass::convert_mass_to_int(weight + (weight / 1000000.0 * upper_mass_limit_ppm as f64)),
            lower_weight_limit: mass::convert_mass_to_int(weight - (weight / 1000000.0 * lower_mass_limit_ppm as f64)),
            thread_count: thread_count,
            decoy_counter: Arc::new(AtomicUsize::new(0)),
            max_modifications_per_decoy: max_modifications_per_decoy,
            fixed_modification_map: Arc::new(fixed_modification_map.clone()),
            variable_modification_map: Arc::new(variable_modification_map.clone()),
            one_amino_acid_substitute_map: Arc::new(*Self::get_one_amino_acid_substitute_map(fixed_modification_map))
        }
    }

    pub fn get_lower_weight_limit(&self) -> i64 {
        return self.lower_weight_limit;
    }

    pub fn get_upper_weight_limit(&self) -> i64 {
        return self.upper_weight_limit;
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

    pub fn generate_decoys(&self, number_of_decoys_to_generate: usize) {
        println!(
            "Need to build {} decoys with weight {} to {}",
            number_of_decoys_to_generate,
            mass::convert_mass_to_float(self.lower_weight_limit),
            mass::convert_mass_to_float(self.upper_weight_limit)
        );
        if number_of_decoys_to_generate > 0 {
            // create threadpoll
            let thread_pool = ThreadPool::new(self.thread_count);
            // loop for starting threads
            let start_time: f64 = time::precise_time_s();
            for _ in 0..self.thread_count {
                // create copies of thread safe pointer of DecoyGenerator which can be move into thread
                let decoy_counter_ptr = self.decoy_counter.clone();
                let fixed_modification_map_ptr = self.fixed_modification_map.clone();
                let variable_modification_map_ptr = self.variable_modification_map.clone();
                // copy primitive attributes of DecoyGenerator which can be moved into thread
                let upper_weight_limit = self.upper_weight_limit;
                let lower_weight_limit = self.lower_weight_limit;
                let max_modifications_per_decoy = self.max_modifications_per_decoy;
                // start thread
                thread_pool.execute(move||{
                    let conn: postgres::Connection = DatabaseConnection::get_database_connection();
                    // create random number generator
                    let mut rng = rand::thread_rng();
                    // get haviest and lightest amino acids
                    let haviest_amino_acid: AminoAcid = AminoAcid::get_haviest();
                    let haviest_amino_acid_weight_with_modification: i64 = match fixed_modification_map_ptr.get(&haviest_amino_acid.get_one_letter_code()){
                        Some(ref modification) => modification.get_mono_mass() + haviest_amino_acid.get_mono_mass(),
                        None => haviest_amino_acid.get_mono_mass()
                    };
                    // let lightest_amino_acid = AminoAcid::get_lightest();
                    // endless loop with label 'decoy_loop
                    'decoy_loop: loop {
                        // create new empty decoy
                        let mut new_decoy: NewDecoy = NewDecoy::new(lower_weight_limit, upper_weight_limit);
                        // get array of fitting amino acids
                        //let distribution_array = *Self::generate_amino_acid_distribution_array();
                        // repeat until new_decoy's weight is smaller than the haviest amino acid.
                        while new_decoy.get_distance_to_mass_tolerance() >= haviest_amino_acid_weight_with_modification {
                            // pick amino acids one letter code at index
                            //let aa_one_letter_code: char = *distribution_array.choose(&mut rng).unwrap();
                            let aa_one_letter_code: char = *AMINO_ACIDS_FOR_DECOY_GENERATION.choose(&mut rng).unwrap();
                            // one letter code to amino acid
                            let random_amino_acid: AminoAcid = AminoAcid::get(aa_one_letter_code);
                            let modification_option = match fixed_modification_map_ptr.get(&random_amino_acid.get_one_letter_code()){
                                Some(ref modification) => Some((*modification).clone()),
                                None => None
                            };
                            new_decoy.push_amino_acid_and_fix_modification(&random_amino_acid, &modification_option);
                            //println!("{} => {}", new_decoy.get_aa_sequence(), mass::convert_mass_to_float(new_decoy.get_weight()));
                        }
                        let mut still_fitting_amino_acids = *Self::get_amino_acids_with_weight_less_or_equals_than(new_decoy.get_distance_to_mass_tolerance(), &fixed_modification_map_ptr);
                        while still_fitting_amino_acids.len() > 0 {
                            // pick amino acids one letter code at index
                            let aa_one_letter_code: char = *still_fitting_amino_acids.choose(&mut rng).unwrap();
                            // one letter code to amino acid
                            let random_amino_acid: AminoAcid = AminoAcid::get(aa_one_letter_code);
                            let modification_option = match fixed_modification_map_ptr.get(&random_amino_acid.get_one_letter_code()){
                                Some(ref modification) => Some((*modification).clone()),
                                None => None
                            };
                            match new_decoy.push_amino_acid_and_fix_modification(&random_amino_acid, &modification_option) {
                                Ok(_) => (),
                                Err(_) => ()
                            }
                            still_fitting_amino_acids = *Self::get_amino_acids_with_weight_less_or_equals_than(new_decoy.get_distance_to_mass_tolerance(), &fixed_modification_map_ptr);
                            //println!("{} => {}", new_decoy.get_aa_sequence(), mass::convert_mass_to_float(new_decoy.get_weight()));
                        }
                        let mut base_decoy: BaseDecoy = new_decoy.as_base_decoy();
                        // check if decoy is peptide
                        match Peptide::exists_where(&conn, "aa_sequence = $1", &[&new_decoy.get_aa_sequence()]) {
                            Ok(_) => continue 'decoy_loop,  // if decoy is peptide, continue with next iteration
                            Err(err) => match err {
                                QueryError::NoMatch => (),  // if no match is found continue with this decoy
                                _ => panic!("proteomic::utility::decoy_generator.generate_decoys() could not check if decoy is peptide: {}", err)
                            }
                        }
                        if new_decoy.hits_mass_tolerance() {
                            match base_decoy.create(&conn) {
                                Ok(query_ok) => match query_ok {
                                    QueryOk::Created => {
                                        let last_count = decoy_counter_ptr.fetch_add(1, Ordering::Relaxed);
                                        println!("count: {}", last_count + 1);
                                    },
                                    QueryOk::AlreadyExists => (),
                                    _ => panic!("proteomic::utility::decoy_generator::generate(): In fact not other QueryOk than QueryOk::Created and QueryOk::AlreadyExists are used in BaseDecoy.create(), so this panic shoud not be reached.")
                                },
                                Err(err) => println!("proteomic::utility::decoy_generator::generate() when second BaseDecoy.create(): {}", err)
                            }
                        }
                        // vary peptide
                    }
                });
            }
            thread_pool.join();
            let stop_time: f64 = time::precise_time_s();
            println!("generate {} decoys in {} s", self.decoy_counter.load(Ordering::Relaxed), stop_time - start_time);
        } else {
            println!("There are plenty of decoys within the given range. nothing to do here.")
        }
    }

    pub fn vary_targets(&self, targets: &Vec<Peptide>) -> Box<HashSet<PlainDecoy>> {
        let conn: postgres::Connection = DatabaseConnection::get_database_connection();
        let mut decoys: HashSet<PlainDecoy> = HashSet::new();
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
                let mut new_decoy: NewDecoy = NewDecoy::from_string(aa_sequence_as_string.as_str(), self.get_lower_weight_limit(), self.get_upper_weight_limit(), self.fixed_modification_map.as_ref());
                // let mut base_decoys = new_decoy.as_base_decoy();
                if new_decoy.hits_mass_tolerance() {
                    // match base_decoys.create(&conn) {
                    //     Ok(query_ok) => match query_ok {
                    //         QueryOk::Created => (),
                    //         QueryOk::AlreadyExists => (),
                    //         _ => ()
                    //     },
                    //     Err(query_err) => panic!("proteomic::utility::decoy_generator::DecoyGenerator.vary_decoy() could not create base decoy from shuffled target: {}", query_err)
                    // }
                    decoys.insert(PlainDecoy::new(new_decoy.get_header().as_str(), new_decoy.get_aa_sequence().as_str(), &new_decoy.get_modified_weight()));
                }
            }
        }
        return Box::new(decoys);
    }

    pub fn vary_new_decoy(&self, new_decoy: &mut NewDecoy) {
        unimplemented!();
    }

    /// calculates a substitution map for swapping amino acids with each other.
    /// so one can lookup which difference in weight a substitution of amino acid x with y has.
    /// this function also considers fixed modifications
    fn get_one_amino_acid_substitute_map(fixed_modification_map: &HashMap<char, Modification>) -> Box<HashMap<char, HashMap<char, i64>>> {
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