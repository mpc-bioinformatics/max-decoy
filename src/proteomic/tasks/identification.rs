use std::collections::{HashMap, HashSet};
use std::io::prelude::*;
use std::io::LineWriter;
use std::fs::OpenOptions;


use proteomic::models::amino_acids::amino_acid::AminoAcid;
use proteomic::models::amino_acids::modification::Modification;
use proteomic::utility::mz_ml::mz_ml_reader::MzMlReader;
use proteomic::utility::decoy_generator::DecoyGenerator;
use proteomic::utility::database_connection::DatabaseConnection;
use proteomic::models::persistable::Persistable;
use proteomic::models::peptides::peptide::Peptide;
use proteomic::models::peptides::peptide_interface::PeptideInterface;
use proteomic::models::peptides::modified_peptide::ModifiedPeptide;
use proteomic::models::peptides::decoy::Decoy;
use proteomic::models::mass;
use proteomic::models::fasta_entry::FastaEntry;


const QUERY_LIMIT_SIZE: i64 = 1000;
const TARGET_DECOY_QUERY_CONDITION: &str = "weight BETWEEN $1 AND $2";


pub struct IdentificationArguments {
    modification_csv_file: String,
    spectrum_file: String,
    max_number_of_variable_modification_per_decoy: u8,
    number_of_decoys: usize,
    lower_mass_tolerance: i64,
    upper_mass_tolerance: i64,
    thread_count: usize
}

impl IdentificationArguments {
    pub fn get_modification_csv_file(&self) -> &str {
        return self.modification_csv_file.as_str();
    }

    pub fn get_spectrum_file(&self) -> &str {
        return self.spectrum_file.as_str();
    }

    pub fn get_max_number_of_variable_modification_per_decoy(&self) -> u8 {
        return self.max_number_of_variable_modification_per_decoy;
    }

    pub fn get_number_of_decoys(&self) -> usize {
        return self.number_of_decoys;
    }

    pub fn get_lower_mass_tolerance(&self) -> i64 {
        return self.lower_mass_tolerance;
    }

    pub fn get_upper_mass_tolerance(&self) -> i64 {
        return self.upper_mass_tolerance;
    }

    pub fn get_thread_count(&self) -> usize {
        return self.thread_count;
    }

    pub fn from_cli_args(cli_args: &clap::ArgMatches) -> IdentificationArguments {
        let modification_csv_file: &str = match cli_args.value_of("MODIFICATION_FILE") {
            Some(modification_csv_file) => modification_csv_file,
            None => panic!("proteomic::tasks::identification::parse_identification_cli_arguments(): you must specify a modification-file")
        };
        let spectrum_file: &str = match cli_args.value_of("SPECTRUM_FILE") {
            Some(spectrum_file) => spectrum_file,
            None => panic!("proteomic::tasks::identification::parse_identification_cli_arguments(): you must specify a spectrum-file")
        };
        let max_number_of_variable_modification_per_decoy: u8 = match cli_args.value_of("MAX_NUMBER_OF_VARIABLE_MODIFICATION_PER_PEPTIDE") {
            Some(number_string) => match number_string.to_owned().parse::<u8>() {
                Ok(number) => number,
                Err(_) => panic!("proteomic::tasks::identification::parse_identification_cli_arguments(): could not cast number-of-decoys to unsigned integer")
            },
            None => 0
        };
        let number_of_decoys: usize = match cli_args.value_of("NUMBER_OF_DECOYS") {
            Some(number_string) => match number_string.to_owned().parse::<usize>() {
                Ok(number) => number,
                Err(_) => panic!("proteomic::tasks::identification::parse_identification_cli_arguments(): could not cast number-of-decoys to unsigned integer")
            },
            None => 1000
        };
        let lower_mass_tolerance: i64 = match cli_args.value_of("LOWER_MASS_TOLERANCE") {
            Some(number_string) => match number_string.to_owned().parse::<i64>() {
                Ok(number) => number,
                Err(_) => panic!("proteomic::tasks::identification::parse_identification_cli_arguments(): could not cast lower-mass-tolerance to integer")
            },
            None => 5
        };
        let upper_mass_tolerance: i64 = match cli_args.value_of("UPPER_MASS_TOLERANCE") {
            Some(number_string) => match number_string.to_owned().parse::<i64>() {
                Ok(number) => number,
                Err(_) => panic!("proteomic::tasks::identification::parse_identification_cli_arguments(): could not cast upper-mass-tolerance to integer")
            },
            None => 5
        };
        let thread_count: usize = match cli_args.value_of("THREAD_COUNT") {
            Some(number_string) => match number_string.to_owned().parse::<usize>() {
                Ok(number) => number,
                Err(_) => panic!("proteomic::tasks::identification::parse_identification_cli_arguments(): could not cast thread-count to unsigned integer")
            },
            None => 2
        };
        return Self {
            modification_csv_file: modification_csv_file.to_owned(),
            spectrum_file: spectrum_file.to_owned(),
            max_number_of_variable_modification_per_decoy: max_number_of_variable_modification_per_decoy,
            number_of_decoys: number_of_decoys,
            lower_mass_tolerance: lower_mass_tolerance,
            upper_mass_tolerance: upper_mass_tolerance,
            thread_count: thread_count
        }
    }
}

pub fn identification_task(identification_args: &IdentificationArguments) {
    let conn = DatabaseConnection::get_database_connection();
    // prepare modifications
    let mods = Modification::create_from_csv_file(identification_args.get_modification_csv_file());
    let mut fixed_modifications_map: HashMap<char, Modification> = HashMap::new();
    let mut variable_modifications_map: HashMap<char, Modification> = HashMap::new();
    for modification in mods.iter() {
        if modification.is_fix() {
            fixed_modifications_map.insert(modification.get_amino_acid_one_letter_code(), modification.clone());
        } else {
            variable_modifications_map.insert(modification.get_amino_acid_one_letter_code(), modification.clone());
        }
    }
    let mut sorted_modifyable_amino_acids: Vec<char> = fixed_modifications_map.keys().map(|key| *key).collect::<Vec<char>>();
    sorted_modifyable_amino_acids.sort();
    // prepare condition for target and decoys
    let mut target_decoy_condition = TARGET_DECOY_QUERY_CONDITION.to_owned();
    if sorted_modifyable_amino_acids.len() > 0 {
        target_decoy_condition.push_str(" AND ");
        let mut conditions: Vec<String> = Vec::new();
        for (idx, amino_acid_one_letter_code) in sorted_modifyable_amino_acids.iter().enumerate() {
            conditions.push(format!("{}_count = ${}", amino_acid_one_letter_code, idx + 3));
        }
        target_decoy_condition.push_str(conditions.join(" AND ").as_str());
    }
    target_decoy_condition.push_str(";");
    // initialize mzML-Reader
    let mz_ml_reader = MzMlReader::new(identification_args.get_spectrum_file());
    let spectra = *mz_ml_reader.get_ms_two_spectra();
    // loop through spectra
    for spectrum in spectra.iter() {
        // calculate tolerances and precursor tolerance
        let tolerances = mass::calculate_upper_and_lower_tolerance(
            *spectrum.get_precurso_mass(),
            identification_args.get_upper_mass_tolerance(),
            identification_args.get_lower_mass_tolerance()
        );
        let precursor_tolerance = mass::calculate_precursor_tolerance(
            *spectrum.get_precurso_mass(),
            identification_args.get_upper_mass_tolerance(),
            identification_args.get_lower_mass_tolerance()
        );
        // build array with maximal
        let mut max_modification_counts: HashMap<char, i16> = HashMap::new();
        for amino_acid_one_letter_code in sorted_modifyable_amino_acids.iter() {
            let amino_acid = AminoAcid::get(*amino_acid_one_letter_code);
            let max_modification_count: i16 = match fixed_modifications_map.get(&amino_acid_one_letter_code) {
                Some(modification) => (*spectrum.get_precurso_mass() / (amino_acid.get_mono_mass() + modification.get_mono_mass())) as i16,
                None => continue
            };
            max_modification_counts.insert(*amino_acid_one_letter_code, max_modification_count);
        }
        // get condition values
        let mut condition_values: Vec<(i64, i64, Vec<i16>)> = Vec::new();
        get_max_modifyable_amino_acid_counts(&fixed_modifications_map, *spectrum.get_precurso_mass(), &tolerances, &sorted_modifyable_amino_acids, &max_modification_counts, 0, &mut Vec::new(), &mut condition_values);
        // gether targets
        println!("search targets...");
        let mut targets: HashSet<FastaEntry> = HashSet::new();
        let mut decoys: HashSet<FastaEntry> = HashSet::new();
        let mut start_time: f64 = time::precise_time_s();
        for query_values in condition_values {
            let mut values: Vec<&postgres::types::ToSql> = Vec::new();
            values.push(&query_values.0);
            values.push(&query_values.1);
            for count in query_values.2.iter() {
                values.push(count);
            }
            let possible_targets = match Peptide::find_where(&conn, target_decoy_condition.as_str(), values.as_ref()) {
                Ok(targets) => targets,
                Err(err) => panic!("proteomic::tasks::identification::identification_task(): could not gether targets: {}", err)
            };
            for peptide in possible_targets {
                #[allow(unused_assignments)] // `modified_target_fits_precursor_tolerance` is actually read in if-instruction below
                let mut modified_target_fits_precursor_tolerance = false;
                let mut modified_peptide = ModifiedPeptide::from_peptide(&peptide, *spectrum.get_precurso_mass(),  precursor_tolerance.0,  precursor_tolerance.1, &fixed_modifications_map);
                modified_target_fits_precursor_tolerance = modified_peptide.hits_mass_tolerance();
                if !modified_target_fits_precursor_tolerance {
                    modified_target_fits_precursor_tolerance = modified_peptide.try_variable_modifications(identification_args.get_max_number_of_variable_modification_per_decoy(), &variable_modifications_map);
                }
                if modified_target_fits_precursor_tolerance {
                    targets.insert(
                        FastaEntry::new(
                            peptide.get_header_with_modification_summary(modified_peptide.get_modification_summary_for_header().as_str()).as_str(),
                            peptide.get_aa_sequence()
                        )
                    );
                }
            }
            if decoys.len() < identification_args.get_number_of_decoys() {
                let mut possible_decoys = match Decoy::find_where(&conn, target_decoy_condition.as_str(), values.as_ref()) {
                    Ok(count) => count,
                    Err(err) => panic!("proteomic::tasks::identification::identification_task(): could not gether decoy: {}", err)
                };
                for decoy in possible_decoys.iter_mut() {
                    #[allow(unused_assignments)] // `modified_decoys_fits_precursor_tolerance` is actually read in if-instruction below
                    let mut modified_decoys_fits_precursor_tolerance = false;
                    let mut modified_decoy = ModifiedPeptide::from_decoy(&decoy, *spectrum.get_precurso_mass(),  precursor_tolerance.0,  precursor_tolerance.1, &fixed_modifications_map);
                    modified_decoys_fits_precursor_tolerance = modified_decoy.hits_mass_tolerance();
                    if !modified_decoys_fits_precursor_tolerance {
                        modified_decoys_fits_precursor_tolerance = modified_decoy.try_variable_modifications(identification_args.get_max_number_of_variable_modification_per_decoy(), &variable_modifications_map);
                    }
                    if modified_decoys_fits_precursor_tolerance {
                        decoy.set_modification_summary(modified_decoy.get_modification_summary_for_header().as_str());
                        decoys.insert(
                            FastaEntry::new(
                                decoy.get_header().as_str(),
                                decoy.get_aa_sequence()
                            )
                        );
                    }
                    if decoys.len() >= identification_args.get_number_of_decoys() { break; }
                }
            }
        }
        let mut stop_time: f64 = time::precise_time_s();
        println!("found {} targets and {} decoys in {} s\nsearching decoys in database...", targets.len(), decoys.len(), stop_time - start_time);
        // gether decoys
        if decoys.len() < identification_args.get_number_of_decoys()  {
            println!("generating decoys...");
            let mut remaining_number_of_decoys = identification_args.get_number_of_decoys() - decoys.len();
            let remaining_number_of_decoys_for_output = remaining_number_of_decoys;
            let generator: DecoyGenerator = DecoyGenerator::new(
                *spectrum.get_precurso_mass(),
                precursor_tolerance.0,
                precursor_tolerance.1,
                identification_args.get_thread_count(),
                identification_args.get_max_number_of_variable_modification_per_decoy(),
                &fixed_modifications_map,
                &variable_modifications_map
            );
            start_time = time::precise_time_s();
            generator.generate_decoys(remaining_number_of_decoys);
            match generator.get_decoys().lock() {
                Ok(generated_decoys) => {
                    for decoy in generated_decoys.iter() {
                        decoys.insert(
                            FastaEntry::new(
                                decoy.get_header().as_str(),
                                decoy.get_aa_sequence()
                            )
                        );
                    }
                }
                Err(_) => panic!("proteomic::tasks::identification::identification_task(): try to lock poisened mutex for decoys at generator.get_decoys().lock()")
            };
            stop_time = time::precise_time_s();
            println!("generate {} decoys in {} s", remaining_number_of_decoys_for_output, stop_time - start_time);
        }
        let fasta_file = match OpenOptions::new().read(true).write(true).create(true).open(format!("{}.fasta", spectrum.get_title()).as_str()) {
            Ok(file) => file,
            Err(err) => panic!("proteomic::tasks::identification::identification_task(): error at opening fasta-file: {}", err)
        };
        let mut fasta_file = LineWriter::new(fasta_file);
        for fasta_entry in targets {
            match fasta_file.write(fasta_entry.to_string().as_bytes()) {
                Ok(_) => (),
                Err(err) => println!("proteomic::tasks::identification::identification_task(): Could not write to file: {}", err)
            }
        }
        for fasta_entry in decoys {
            match fasta_file.write(fasta_entry.to_string().as_bytes()) {
                Ok(_) => (),
                Err(err) => println!("proteomic::tasks::identification::identification_task(): Could not write to file: {}", err)
            }
        }
    }
}

fn get_max_modifyable_amino_acid_counts(fixed_modifications_map: &HashMap<char, Modification>, precursor: i64, tolerances: &(i64, i64), amino_acid_one_letter_codes: &Vec<char>, max_modification_counts: &HashMap<char, i16>, amino_acid_index: usize, count_combination: &mut Vec<i16>, results: &mut Vec<(i64, i64, Vec<i16>)>) {
    let amino_acid_one_letter_code = match amino_acid_one_letter_codes.get(amino_acid_index) {
        Some(one_letter_code) => one_letter_code,
        None => &'_'
    };
    if let Some(max_modification_count) = max_modification_counts.get(amino_acid_one_letter_code) {
        if let Some(ref modification) = fixed_modifications_map.get(amino_acid_one_letter_code) {
            let mut new_precursor = precursor;
            for mod_count in 0..*max_modification_count {
                new_precursor -= modification.get_mono_mass();
                if precursor > 0 {
                    count_combination.push(mod_count);
                    if amino_acid_index < max_modification_counts.len() - 1 {
                        get_max_modifyable_amino_acid_counts(fixed_modifications_map, new_precursor, tolerances,  amino_acid_one_letter_codes, max_modification_counts, amino_acid_index + 1, count_combination, results);
                    } else {
                        results.push((
                            precursor - tolerances.0,
                            precursor + tolerances.1,
                            count_combination.clone()
                        ));
                    }
                    count_combination.pop();
                }
            }
        }
    }
}