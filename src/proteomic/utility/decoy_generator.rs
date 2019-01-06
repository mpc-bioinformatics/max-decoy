extern crate rand;
extern crate threadpool;


use std::collections::HashMap;
use std::sync::Mutex;
use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering};

use self::rand::seq::SliceRandom;
use self::threadpool::ThreadPool;

use proteomic::models::mass;
use proteomic::utility::database_connection::DatabaseConnection;
use proteomic::utility::combinations::n_choose_k::NChooseK;
// use proteomic::utility::mass::NeutralLoss;
// use proteomic::models::peptide::Peptide;
use proteomic::models::persistable::{Persistable, QueryOk};
use proteomic::models::decoys::decoy::{Decoy, PlainDecoy};
use proteomic::models::decoys::new_decoy::{NewDecoy, NewDecoyError};
use proteomic::models::decoys::base_decoy::BaseDecoy;
use proteomic::models::amino_acids::amino_acid::AminoAcid;
use proteomic::models::amino_acids::modification::Modification;

const AMINO_ACIDS_FOR_DECOY_GENERATION: &'static [char] = &['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'J', 'K', 'M', 'F', 'P', 'O', 'S', 'T', 'V', 'W', 'Y'];

pub struct DecoyGenerator {
    upper_weight_limit: i64,
    lower_weight_limit: i64,
    number_of_decoys_to_generate: usize,
    thread_count: usize,
    decoy_counter: Arc<AtomicUsize>,
    max_modifications_per_decoy: i32,
    fixed_modification_map: Arc<Mutex<HashMap<char, Modification>>>,
    variable_modification_map: Arc<Mutex<HashMap<char, Modification>>>
}

impl DecoyGenerator {
    pub fn new(weight: f64, lower_mass_limit_ppm: i64, upper_mass_limit_ppm: i64, number_of_decoys_to_generate: usize, thread_count: usize, max_modifications_per_decoy: i32, fixed_modification_map: &HashMap<char, Modification>, variable_modification_map: &HashMap<char, Modification>) -> Self {
        let upper_weight_limit = mass::convert_mass_to_int(weight + (weight / 1000000.0 * upper_mass_limit_ppm as f64));
        let lower_weight_limit = mass::convert_mass_to_int(weight - (weight / 1000000.0 * lower_mass_limit_ppm as f64));
        let conn = DatabaseConnection::get_database_connection();
        let already_fitting_decoys: usize = match PlainDecoy::count_where_mass_tolerance(&conn, lower_weight_limit, upper_weight_limit) {
            Ok(count) => count as usize,
            Err(_) => 0
        };
        // make sure that only the difference between number of decoys in database and number of requested decoys are generated
        let mut decoy_difference: usize = 0;
        if already_fitting_decoys < number_of_decoys_to_generate {
            decoy_difference = number_of_decoys_to_generate - already_fitting_decoys;
        }
        println!(
            "Should generate {} decoys within mass tolerance of {} to {}.",
            number_of_decoys_to_generate,
            mass::convert_mass_to_float(lower_weight_limit), 
            mass::convert_mass_to_float(upper_weight_limit)
        );
        println!("Found {} already fitting decoys, will generate {} more.", already_fitting_decoys, decoy_difference);
        return DecoyGenerator{
            upper_weight_limit: upper_weight_limit,
            lower_weight_limit: lower_weight_limit,
            number_of_decoys_to_generate: decoy_difference,
            thread_count: thread_count,
            decoy_counter: Arc::new(AtomicUsize::new(0)),
            max_modifications_per_decoy: max_modifications_per_decoy,
            fixed_modification_map: Arc::new(Mutex::new(fixed_modification_map.clone())),
            variable_modification_map: Arc::new(Mutex::new(variable_modification_map.clone()))
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
    fn get_amino_acids_with_weight_less_or_equals_than(weight: i64) -> Box<Vec<char>> {
        let mut amino_acids: Vec<char> = Vec::new();
        for aa_one_letter_code in AMINO_ACIDS_FOR_DECOY_GENERATION {
            let current_amino_acid: AminoAcid = AminoAcid::get(*aa_one_letter_code);
            if current_amino_acid.get_mono_mass() <= weight {
                amino_acids.push(current_amino_acid.get_one_letter_code());
            }
        }
        return Box::new(amino_acids);
    }

    pub fn generate_decoys(&self) {
        if self.number_of_decoys_to_generate > 0 {
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
                let number_of_decoys_to_generate = self.number_of_decoys_to_generate;
                let max_modifications_per_decoy = self.max_modifications_per_decoy;
                // start thread
                thread_pool.execute(move||{
                    let conn: postgres::Connection = DatabaseConnection::get_database_connection();
                    // create random number generator
                    let mut rng = rand::thread_rng();
                    // get haviest and lightest amino acids
                    let haviest_amino_acid: AminoAcid = AminoAcid::get_haviest();
                    // let lightest_amino_acid = AminoAcid::get_lightest();
                    // endless loop with label 'decoy_loop
                    'decoy_loop: loop {
                        // create new empty decoy
                        let mut new_decoy: NewDecoy = NewDecoy::new(lower_weight_limit, upper_weight_limit);
                        // get array of fitting amino acids
                        let distribution_array = *Self::generate_amino_acid_distribution_array();
                        // repeat until new_decoy's weight is smaller than the haviest amino acid.
                        while new_decoy.get_distance_to_hit_mass_tolerance() >= haviest_amino_acid.get_mono_mass() {
                            // pick amino acids one letter code at index
                            let aa_one_letter_code: char = *distribution_array.choose(&mut rng).unwrap();
                            // one letter code to amino acid
                            let random_amino_acid: AminoAcid = AminoAcid::get(aa_one_letter_code);
                            new_decoy.push_amino_acid(&random_amino_acid);
                            //println!("{} => {}", new_decoy.get_aa_sequence(), new_decoy.get_weight());
                        }
                        let mut still_fitting_amino_acids = *Self::get_amino_acids_with_weight_less_or_equals_than(new_decoy.get_distance_to_hit_mass_tolerance());
                        while still_fitting_amino_acids.len() > 0 {
                            // pick amino acids one letter code at index
                            let aa_one_letter_code: char = *still_fitting_amino_acids.choose(&mut rng).unwrap();
                            // one letter code to amino acid
                            let random_amino_acid: AminoAcid = AminoAcid::get(aa_one_letter_code);
                            new_decoy.push_amino_acid(&random_amino_acid);
                            //println!("{} => {}", new_decoy.get_aa_sequence(), new_decoy.get_weight());
                            still_fitting_amino_acids = *Self::get_amino_acids_with_weight_less_or_equals_than(new_decoy.get_distance_to_hit_mass_tolerance());
                        }
                        let mut base_decoy: BaseDecoy = new_decoy.as_base_decoy();
                        if new_decoy.hits_mass_tolerance() {
                            match base_decoy.create(&conn) {
                                Ok(query_ok) => match query_ok {
                                    QueryOk::Created => {
                                        decoy_counter_ptr.fetch_add(1, Ordering::Relaxed);
                                    },
                                    QueryOk::AlreadyExists => (),
                                    _ => panic!("proteomic::utility::decoy_generator::generate(): In fact not other QueryOk than QueryOk::Created and QueryOk::AlreadyExists are used in BaseDecoy.create(), so this panic shoud not be reached.")
                                },
                                Err(err) => println!("proteomic::utility::decoy_generator::generate() when second BaseDecoy.create(): {}", err)
                            }
                        }
                        if max_modifications_per_decoy > 0 {
                            let mut modification_positions: Vec<i32> = Vec::new();
                            modification_positions.push(-1); // -1 is n_terminus
                            for idx in 0..new_decoy.get_length() {
                                modification_positions.push(idx);
                            }
                            modification_positions.push(-2); // -2 is c_terminus
                            for current_modification_max in 1..max_modifications_per_decoy + 1 {
                                let combinations = NChooseK::new(current_modification_max, modification_positions.clone());
                                for combination in combinations {
                                    for modification_position in combination {
                                        let amino_acid_one_letter_code = match modification_position {
                                            -1 => new_decoy.get_n_terminus_amino_acid(),
                                            -2 => new_decoy.get_c_terminus_amino_acid(),
                                            _ => new_decoy.get_amino_acid_at(modification_position as usize)
                                        };
                                        let mut modification_option: Option<Modification> = match fixed_modification_map_ptr.lock() {
                                            Ok(fixed_modifications_map) => {
                                                match fixed_modifications_map.get(&amino_acid_one_letter_code) {
                                                    Some(modification) => Some(modification.clone()),
                                                    None => None
                                                }
                                            },
                                            Err(_) => panic!("ERROR [DecoyGenerator]: could not gather lock for fixed modifications")
                                        };
                                        if modification_option.is_none() {
                                            modification_option = match variable_modification_map_ptr.lock() {
                                                Ok(variable_modifications_map) => {
                                                    match variable_modifications_map.get(&amino_acid_one_letter_code){
                                                        Some(modification) => Some(modification.clone()),
                                                        None => None
                                                    }
                                                },
                                                Err(_) => panic!("ERROR [DecoyGenerator]: could not gather lock for variable modifications")
                                            };
                                        }
                                        if let Some(modification) = modification_option {
                                            match modification_position {
                                                -1 => match new_decoy.set_n_terminus_modification(&modification) {
                                                    Ok(_) => (),
                                                    Err(new_decoy_err) => match new_decoy_err {
                                                        NewDecoyError::ModificationDoesNotMatchToAminoAcid => panic!("proteomic::utility::decoy_generator::generate(): Amino acid with non-matching modification is passed to NewDecoy.set_n_terminus_modification(). This should not happen here, because the code get matching modification from HashMap with amino acid one letter code.")
                                                    }
                                                },
                                                -2 => match new_decoy.set_c_terminus_modification(&modification) {
                                                    Ok(_) => (),
                                                    Err(new_decoy_err) => match new_decoy_err {
                                                        NewDecoyError::ModificationDoesNotMatchToAminoAcid => panic!("proteomic::utility::decoy_generator::generate(): Amino acid with non-matching modification is passed to NewDecoy.set_c_terminus_modification(). This should not happen here, because the code get matching modification from HashMap with amino acid one letter code.")
                                                    }
                                                },
                                                _ => match new_decoy.set_modification_at(modification_position as usize, &modification) {
                                                    Ok(_) => (),
                                                    Err(new_decoy_err) => match new_decoy_err {
                                                        NewDecoyError::ModificationDoesNotMatchToAminoAcid => panic!("proteomic::utility::decoy_generator::generate(): Amino acid with non-matching modification is passed to NewDecoy.set_modification_at(). This should not happen here, because the code get matching modification from HashMap with amino acid one letter code.")
                                                    }
                                                }
                                            };
                                        }
                                    }
                                    if (new_decoy.get_number_of_modifications() > 0) & new_decoy.hits_mass_tolerance() {
                                        if !base_decoy.is_persisted() {
                                            match base_decoy.create(&conn) {
                                                Ok(query_ok) => match query_ok {
                                                    QueryOk::Created => (),             // do not increase counter here, because BaseDecoy is only created for reference to ModifiedDecoy
                                                    QueryOk::AlreadyExists => (),
                                                    _ => panic!("proteomic::utility::decoy_generator::generate(): In fact not other QueryOk than QueryOk::Created and QueryOk::AlreadyExists are used in BaseDecoy.create(), so this panic shoud not be reached.")
                                                },
                                                Err(err) => println!("proteomic::utility::decoy_generator::generate() when second BaseDecoy.create(): {}", err)
                                            }
                                        }
                                        let mut modified_decoy = new_decoy.as_modified_decoy(&base_decoy);
                                        match modified_decoy.create(&conn) {
                                            Ok(query_ok) => match query_ok {
                                                QueryOk::Created => {
                                                    decoy_counter_ptr.fetch_add(1, Ordering::Relaxed);
                                                },
                                                QueryOk::AlreadyExists => (),
                                                _ => panic!("proteomic::utility::decoy_generator::generate(): In fact not other QueryOk than QueryOk::Created and QueryOk::AlreadyExists are used in ModifiedDecoy.create(), so this panic shoud not be reached.")
                                            },
                                            Err(err) => println!("proteomic::utility::decoy_generator::generate() when ModifiedDecoy.create(): {}", err)
                                        }
                                    }
                                    new_decoy.remove_all_modifications();
                                }
                            }
                            if decoy_counter_ptr.load(Ordering::Relaxed) >= number_of_decoys_to_generate {
                                break;
                            }
                        }
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
}