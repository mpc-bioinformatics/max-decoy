use std::collections::HashMap;
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


pub struct IdentificationArguments {
    modification_csv_file: String,
    spectrum_file: String,
    max_number_of_variable_modification_per_decoy: u8,
    number_of_decoys_per_target: usize,
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

    pub fn get_number_of_decoys_per_target(&self) -> usize {
        return self.number_of_decoys_per_target;
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
        let number_of_decoys_per_target: usize = match cli_args.value_of("NUMBER_OF_DECOYS_PER_TARGET") {
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
            number_of_decoys_per_target: number_of_decoys_per_target,
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
    let mut haviest_fixed_modification: Option<&Modification> = None;
    for (_, modification) in fixed_modifications_map.iter() {
        haviest_fixed_modification = match haviest_fixed_modification {
            Some(ref mut haviest_modification) if haviest_modification.get_mono_mass() < modification.get_mono_mass() => Some(modification),
            Some(_) => continue,
            None => Some(modification)
        };
    }
    let mz_ml_reader = MzMlReader::new(identification_args.get_spectrum_file());
    let spectra = *mz_ml_reader.get_ms_two_spectra();
    for spectrum in spectra.iter() {
        let generator: DecoyGenerator = DecoyGenerator::new(
            *spectrum.get_precurso_mass(),
            identification_args.get_upper_mass_tolerance(),
            identification_args.get_lower_mass_tolerance(),
            identification_args.get_thread_count(),
            identification_args.get_max_number_of_variable_modification_per_decoy(),
            &fixed_modifications_map,
            &variable_modifications_map
        );
        let highest_possiple_sequence_length = *spectrum.get_precurso_mass() / AminoAcid::get_lightest().get_mono_mass() + 1;
        let lower_mass_tolerance_without_fixed_modifications = match haviest_fixed_modification {
            Some(modification) => generator.get_lower_weight_limit() - (highest_possiple_sequence_length * modification.get_mono_mass()),
            None => generator.get_lower_weight_limit()
        };
        let fasta_file = match OpenOptions::new().read(true).write(true).create(true).open(format!("{}.fasta", spectrum.get_title()).as_str()) {
            Ok(file) => file,
            Err(err) => panic!("proteomic::tasks::identification::identification_task(): error at opening fasta-file: {}", err)
        };
        let mut fasta_file = LineWriter::new(fasta_file);
        let mut offset_counter = 0;
        let mut target_counter = 0;
        // gether targets
        loop {
            let possible_targets = match Peptide::find_where(&conn, "weight BETWEEN $1 AND $2 OFFSET $3 LIMIT $4", &[&lower_mass_tolerance_without_fixed_modifications, &generator.get_upper_weight_limit(), &(offset_counter * 1000), &1000]) {
                Ok(targets) => targets,
                Err(err) => panic!("proteomic::tasks::identification::identification_task(): could not gether targets: {}", err)
            };
            if possible_targets.len() == 0 { break; }
            for peptide in possible_targets {
                let mut modified_decoys_fits_precursor_tolerance = false;
                let mut modified_peptide = ModifiedPeptide::from_peptide(&peptide, *spectrum.get_precurso_mass(),  generator.get_lower_weight_limit(),  generator.get_upper_weight_limit(), &fixed_modifications_map);
                if modified_peptide.hits_mass_tolerance() {
                    modified_decoys_fits_precursor_tolerance = true;
                } else {
                    modified_peptide.try_variable_modifications(identification_args.get_max_number_of_variable_modification_per_decoy(), &variable_modifications_map);
                    if modified_peptide.hits_mass_tolerance() { modified_decoys_fits_precursor_tolerance = true; }
                }
                if modified_decoys_fits_precursor_tolerance {
                    match fasta_file.write(
                        Peptide::as_fasta_entry(
                            peptide.get_header_with_modification_summary(&conn, modified_peptide.get_modification_summary_for_header().as_str()).as_str(),
                            peptide.get_aa_sequence()
                        ).as_bytes()
                    ) {
                        Ok(_) => target_counter += 1,
                        Err(err) => println!("proteomic::tasks::identification::identification_task(): Could not write to file: {}", err)
                    }
                }
            }
            offset_counter += 1;
        }
        // gether decoys
        let mut number_of_decyos_to_generate = identification_args.get_number_of_decoys_per_target() * target_counter;
        offset_counter = 0;
        loop {
            let mut possible_decoys = match Decoy::find_where(&conn, "weight BETWEEN $1 AND $2 OFFSET $3 LIMIT $4", &[&lower_mass_tolerance_without_fixed_modifications, &generator.get_upper_weight_limit(), &(offset_counter * 1000), &1000]) {
                Ok(count) => count,
                Err(err) => panic!("proteomic::tasks::identification::identification_task(): could not gether decoy: {}", err)
            };
            if possible_decoys.len() == 0 { break; }
            for decoy in possible_decoys.iter_mut() {
                let mut modified_decoys_fits_precursor_tolerance = false;
                let mut modified_decoy = ModifiedPeptide::from_decoy(&decoy, *spectrum.get_precurso_mass(),  generator.get_lower_weight_limit(),  generator.get_upper_weight_limit(), &fixed_modifications_map);
                modified_decoys_fits_precursor_tolerance = modified_decoy.hits_mass_tolerance();
                if modified_decoys_fits_precursor_tolerance {
                    modified_decoys_fits_precursor_tolerance = modified_decoy.try_variable_modifications(identification_args.get_max_number_of_variable_modification_per_decoy(), &variable_modifications_map);
                }
                if modified_decoys_fits_precursor_tolerance {
                    decoy.set_modification_summary(modified_decoy.get_modification_summary_for_header().as_str());
                    match fasta_file.write(
                        Decoy::as_fasta_entry(
                            decoy.get_header().as_str(),
                            decoy.get_aa_sequence()
                        ).as_bytes()
                    ) {
                        Ok(_) => number_of_decyos_to_generate -= 1,
                        Err(err) => println!("proteomic::tasks::identification::identification_task(): Could not write to file: {}", err)
                    }
                }
            }
            if number_of_decyos_to_generate == 0 { break; }
            offset_counter += 1;
        }
        generator.generate_decoys(number_of_decyos_to_generate);
        match generator.get_decoys().lock() {
            Ok(decoys) => {
                for decoy in decoys.iter() {
                    match fasta_file.write(
                        Decoy::as_fasta_entry(
                            decoy.get_header().as_str(),
                            decoy.get_aa_sequence()
                        ).as_bytes()
                    ) {
                        Ok(_) => number_of_decyos_to_generate -= 1,
                        Err(err) => println!("proteomic::tasks::identification::identification_task(): Could not write to file: {}", err)
                    }
                }
            }
            Err(_) => panic!("proteomic::tasks::identification::identification_task(): try to lock poisened mutex for decoys at generator.get_decoys().lock()")
        };
    }
}