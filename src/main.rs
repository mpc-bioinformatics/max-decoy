extern crate clap;
extern crate num_cpus;
extern crate postgres;
extern crate postgres_array;
extern crate quick_xml;
extern crate sha1;
extern crate rand;
extern crate threadpool;
extern crate onig;
extern crate time;
extern crate csv;
extern crate dotenv;

use std::collections::HashMap;

use clap::{Arg, App, SubCommand};

mod proteomic;
use proteomic::utility::decoy_generator::DecoyGenerator;
use proteomic::utility::mz_ml::mz_ml_reader::MzMlReader;
use proteomic::utility::mz_ml::spectrum::Spectrum;

use proteomic::tasks::indentification::{identification_task, IdentificationArguments};
use proteomic::tasks::digestion::{digest_to_database_task, DigestionArguments};

use proteomic::models::amino_acids::modification::Modification;
use proteomic::models::mass;


fn run_decoy_generation(decoy_generation_cli_args: &clap::ArgMatches) {
    let modification_csv_file: String = match decoy_generation_cli_args.value_of("MODIFICATION_FILE") {
        Some(modification_csv_file) => modification_csv_file.to_owned(),
        None => String::new()
    };
    let max_modifications_per_decoy: u8 = match decoy_generation_cli_args.value_of("MAX_MODIFICATION_PER_DECOY") {
        Some(number_string) => {
            match number_string.to_owned().parse::<u8>() {
                Ok(number) => number,
                Err(_) => panic!("ERROR [decoy-generation]: could not cast max-modification-per-decoy to (unsigned) integer (8 bit)")
            }
        },
        None => 0
    };
    let number_of_decoys: usize = match decoy_generation_cli_args.value_of("NUMBER_OF_DECOYS") {
        Some(number_string) => {
            match number_string.to_owned().parse::<usize>() {
                Ok(number) => number,
                Err(_) => panic!("ERROR [decoy-generation]: could not cast number-of-decoys to unsigned integer")
            }
        },
        None => 1000
    };
    let lower_mass_tolerance: i64 = match decoy_generation_cli_args.value_of("LOWER_MASS_TOLERANCE") {
        Some(number_string) => {
            match number_string.to_owned().parse::<i64>() {
                Ok(number) => number,
                Err(_) => panic!("ERROR [decoy-generation]: could not cast lower-mass-tolerance to integer")
            }
        },
        None => 5
    };
    let upper_mass_tolerance: i64 = match decoy_generation_cli_args.value_of("UPPER_MASS_TOLERANCE") {
        Some(number_string) => {
            match number_string.to_owned().parse::<i64>() {
                Ok(number) => number,
                Err(_) => panic!("ERROR [decoy-generation]: could not cast upper-mass-tolerance to integer")
            }
        },
        None => 5
    };
    let precursor_mass: i64 = match decoy_generation_cli_args.value_of("PRECURSOR_MASS") {
        Some(number_string) => {
            match number_string.to_owned().parse::<f64>() {
                Ok(number) => mass::convert_mass_to_int(number),
                Err(_) => panic!("ERROR [decoy-generation]: could not cast upper-mass-tolerance to integer")
            }
        },
        None => 5
    };
    let cpu_thread_count: usize = num_cpus::get();
    let thread_count: usize = match decoy_generation_cli_args.value_of("THREAD_COUNT") {
        Some(count) => {
            match count.to_lowercase().as_str() {
                "max" => cpu_thread_count,
                _ => match count.to_owned().parse::<usize>() {
                    Ok(mut parsed_count) => {
                        let cpu_thread_count: usize = num_cpus::get();
                        if (1 > parsed_count) | (parsed_count > num_cpus::get()) {
                            panic!("ERROR [decoy-generation]: threadcount must be between {} and {}.", 1, cpu_thread_count);
                        }
                        parsed_count
                    } ,
                    Err(_err) => {
                        panic!("ERROR [decoy-generation]: cannot parse thread-count to (unsigned) integer");
                    }
                }
            }
        },
        None => {
            println!("WARNING [decoy-generation]: no thread-count , set it to 1.");
            1
        }
    };
    // prepare modifications
    let mods = Modification::create_from_csv_file(&modification_csv_file);
    let mut fixed_modifications_map: HashMap<char, Modification> = HashMap::new();
    let mut variable_modifications_map: HashMap<char, Modification> = HashMap::new();
    for modification in mods.iter() {
        if modification.is_fix() {
            fixed_modifications_map.insert(modification.get_amino_acid_one_letter_code(), modification.clone());
        } else {
            variable_modifications_map.insert(modification.get_amino_acid_one_letter_code(), modification.clone());
        }
    }
    let start_time: f64 = time::precise_time_s();
    let generator: DecoyGenerator = DecoyGenerator::new(precursor_mass, upper_mass_tolerance, lower_mass_tolerance, thread_count, max_modifications_per_decoy, &fixed_modifications_map, &variable_modifications_map);
    generator.generate_decoys(number_of_decoys);
    let stop_time: f64 = time::precise_time_s();
    println!("generate {} decoys in {} s", number_of_decoys, stop_time - start_time);
}


fn run_amino_acid_substitution(substitution_cli_args: &clap::ArgMatches) {
    let modification_csv_file: String = match substitution_cli_args.value_of("MODIFICATION_FILE") {
        Some(modification_csv_file) => modification_csv_file.to_owned(),
        None => String::new()
    };
    let source_amino_acid: char = match substitution_cli_args.value_of("SOURCE_AMINO_ACID") {
        Some(string) => {
            match string.chars().next() {
                Some(character) => character,
                None => panic!("ERROR [decoy-generation]: could not parse source-amino-acid to character")
            }
        },
        None => '_'
    };
    let destination_amino_acid: char = match substitution_cli_args.value_of("DESTINATION_AMINO_ACID") {
        Some(string) => {
            match string.chars().next() {
                Some(character) => character,
                None => panic!("ERROR [decoy-generation]: could not parse destination-amino-acid to character")
            }
        },
        None => '_'
    };
    // prepare modifications
    let mods = Modification::create_from_csv_file(&modification_csv_file);
    let mut fixed_modifications_map: HashMap<char, Modification> = HashMap::new();
    let mut variable_modifications_map: HashMap<char, Modification> = HashMap::new();
    for modification in mods.iter() {
        if modification.is_fix() {
            fixed_modifications_map.insert(modification.get_amino_acid_one_letter_code(), modification.clone());
        } else {
            variable_modifications_map.insert(modification.get_amino_acid_one_letter_code(), modification.clone());
        }
    }
    let substitution_map = *DecoyGenerator::get_one_amino_acid_substitute_map(&fixed_modifications_map);
    if let Some(swaps) = substitution_map.get(&source_amino_acid) {
        if let Some(weight_change) = swaps.get(&destination_amino_acid) {
            println!("{} => {} = {}", source_amino_acid, destination_amino_acid, mass::convert_mass_to_float(*weight_change));
        }
    }
}

/// Splits up the given mzML-file into mzML-files containing a single MS2-spectrum and writes them into the given destination folder
fn run_spectrum_splitup(mz_ml_splitup_cli_args: &clap::ArgMatches) {
    // paesing cli arguments
    let mz_ml_file: &str = match mz_ml_splitup_cli_args.value_of("MZ_ML_FILE") {
        Some(mz_ml_file) => mz_ml_file,
        None => panic!("you must provide a MzMl-file")
    };
    let destination_folder: &str = match mz_ml_splitup_cli_args.value_of("DESTINATION_FOLDER") {
        Some(destination_folder) => destination_folder,
        None => panic!("you must provide a destination folder")
    };
    let file_suffix: &str = match mz_ml_splitup_cli_args.value_of("FILE_SUFFIX") {
        Some(file_suffix) => file_suffix,
        None => ""
    };
    let start_time: f64 = time::precise_time_s();
    let mz_ml_reader = MzMlReader::new(mz_ml_file);
    let content_before_spectrum_list = mz_ml_reader.get_content_before_spectrum_list();
    let spectra: Vec<Spectrum> = *mz_ml_reader.get_ms_two_spectra();
    for spectrum in spectra.iter() {
        spectrum.to_mz_ml(content_before_spectrum_list.as_str(), destination_folder, file_suffix);
    }
    let stop_time: f64 = time::precise_time_s();
    println!("found and write {} MS2-spectra in {} s", spectra.len(), stop_time - start_time);
}

/// Splits up the given mzML-file into mzML-files containing a single MS2-spectrum and writes them into the given destination folder
fn run_identification(cli_args: &clap::ArgMatches) {

}

fn main() {
    let matches = App::new("Peptide Magic")
    .version("1.0")
    .author("Dirk Winkelhardt <dirk.winkelhardt@gmail.com>")
    .about("Does awesome things")
    .subcommand(
        SubCommand::with_name("digest")
        .arg(
            Arg::with_name("INPUT_FILE")
            .short("i")
            .long("input-file")
            .value_name("INPUT_FILE")
            .takes_value(true)
        )
        .arg(
            Arg::with_name("INPUT_FORMAT")
            .short("f")
            .long("format")
            .value_name("FORMAT")
            .takes_value(true)
        )
        .arg(
            Arg::with_name("THREAD_COUNT")
            .short("t")
            .long("thread-count")
            .value_name("THREAD_COUNT")
            .takes_value(true)
        )
        .arg(
            Arg::with_name("NUMBER_OF_MISSED_CLEAVAGES")
            .short("c")
            .long("number-of-missed-cleavages")
            .value_name("NUMBER_OF_MISSED_CLEAVAGES")
            .takes_value(true)
            .help("default: 2, maximal: 60")
        )
        .arg(
            Arg::with_name("MIN_PEPTIDE_LENGTH")
            .short("l")
            .long("minimum-peptide_length")
            .value_name("MIN_PEPTIDE_LENGTH")
            .takes_value(true)
        )
        .arg(
            Arg::with_name("MAX_PEPTIDE_LENGTH")
            .short("h")
            .long("maximum-peptide_length")
            .value_name("MAX_PEPTIDE_LENGTH")
            .takes_value(true)
            .help("maximal: 60")
        )
        .arg(
            Arg::with_name("ENZYM_NAME")
            .short("e")
            .long("enzym-name")
            .value_name("ENZYM_NAME")
            .takes_value(true)
            .help("Trypsin")
        )
    )
    .subcommand(
        SubCommand::with_name("decoy-generation")
        .arg(
            Arg::with_name("MODIFICATION_FILE")
            .short("m")
            .long("modification-file")
            .value_name("INPUT_FILE")
            .takes_value(true)
        )
        .arg(
            Arg::with_name("MAX_MODIFICATION_PER_DECOY")
            .short("n")
            .long("max-modification-per-decoy")
            .value_name("MAX_MODIFICATION_PER_DECOY")
            .takes_value(true)
            .help("Integer, Default: 0")
        )
        .arg(
            Arg::with_name("PRECURSOR_MASS")
            .short("p")
            .long("precursor-mass")
            .value_name("PRECURSOR_MASS")
            .takes_value(true)
        )
        .arg(
            Arg::with_name("NUMBER_OF_DECOYS")
            .short("d")
            .long("number-of-decoys")
            .value_name("NUMBER_OF_DECOYS")
            .takes_value(true)
            .help("Integer, Default: 1000")
        )
        .arg(
            Arg::with_name("LOWER_MASS_TOLERANCE")
            .short("l")
            .long("lower-mass-tolerance")
            .value_name("LOWER_MASS_TOLERANCE")
            .takes_value(true)
            .help("Integer, Unit: ppm, Default: 5")
        )
        .arg(
            Arg::with_name("UPPER_MASS_TOLERANCE")
            .short("u")
            .long("upper-mass-tolerance")
            .value_name("UPPER_MASS_TOLERANCE")
            .takes_value(true)
            .help("Integer, Unit: ppm, Default: 5")
        )
        .arg(
            Arg::with_name("THREAD_COUNT")
            .short("t")
            .long("thread-count")
            .value_name("THREAD_COUNT")
            .takes_value(true)
        )
    )
    .subcommand(
        SubCommand::with_name("spectrum-splitup")
        .arg(
            Arg::with_name("MZ_ML_FILE")
            .short("m")
            .long("mz-ml-file")
            .value_name("MZ_ML_FILE")
            .takes_value(true)
        )
        .arg(
            Arg::with_name("DESTINATION_FOLDER")
            .short("d")
            .long("destination-folder")
            .value_name("DESTINATION_FOLDER")
            .takes_value(true)
        )
        .arg(
            Arg::with_name("FILE_SUFFIX")
            .short("s")
            .long("file-suffix")
            .value_name("FILE_SUFFIX")
            .takes_value(true)
            .help("Helping identifying your mzML among others.")
        )
    )
    .subcommand(
        SubCommand::with_name("amino-acid-substitution")
        .arg(
            Arg::with_name("MODIFICATION_FILE")
            .short("m")
            .long("modification_file")
            .value_name("MODIFICATION_FILE")
            .takes_value(true)
        )
        .arg(
            Arg::with_name("SOURCE_AMINO_ACID")
            .short("s")
            .long("source-amino-acid")
            .value_name("SOURCE_AMINO_ACID")
            .takes_value(true)
        )
        .arg(
            Arg::with_name("DESTINATION_AMINO_ACID")
            .short("d")
            .long("destination-amino-acid")
            .value_name("DESTINATION_AMINO_ACID")
            .takes_value(true)
        )
    )
    .subcommand(
        SubCommand::with_name("identification")
        .arg(
            Arg::with_name("MODIFICATION_FILE")
            .short("m")
            .long("modification-file")
            .value_name("INPUT_FILE")
            .takes_value(true)
        )
        .arg(
            Arg::with_name("SPECTRUM_FILE")
            .short("m")
            .long("spectrum-file")
            .value_name("SPECTRUM_FILE")
            .takes_value(true)
            .help("mzML-file")
        )
        .arg(
            Arg::with_name("MAX_NUMBER_OF_VARIABLE_MODIFICATION_PER_PEPTIDE")
            .short("n")
            .long("max-number-of-variable-modification-per-peptide")
            .value_name("MAX_NUMBER_OF_VARIABLE_MODIFICATION_PER_PEPTIDE")
            .takes_value(true)
            .help("Integer, Default: 0")
        )
        .arg(
            Arg::with_name("NUMBER_OF_DECOYS_PER_TARGET")
            .short("d")
            .long("number-of-decoys-per-target")
            .value_name("NUMBER_OF_DECOYS_PER_TARGET")
            .takes_value(true)
            .help("Integer, Default: 1000")
        )
        .arg(
            Arg::with_name("LOWER_MASS_TOLERANCE")
            .short("l")
            .long("lower-mass-tolerance")
            .value_name("LOWER_MASS_TOLERANCE")
            .takes_value(true)
            .help("Integer, Unit: ppm, Default: 5")
        )
        .arg(
            Arg::with_name("UPPER_MASS_TOLERANCE")
            .short("u")
            .long("upper-mass-tolerance")
            .value_name("UPPER_MASS_TOLERANCE")
            .takes_value(true)
            .help("Integer, Unit: ppm, Default: 5")
        )
        .arg(
            Arg::with_name("THREAD_COUNT")
            .short("t")
            .long("thread-count")
            .value_name("THREAD_COUNT")
            .takes_value(true)
            .help("Integer, Default: 2")
        )
    )
    .get_matches();


    if let Some(cli_args) = matches.subcommand_matches("digest") {
        let digestion_args = DigestionArguments::from_cli_args(cli_args);
        digest_to_database_task(&digestion_args);
    }
    if let Some(cli_args) = matches.subcommand_matches("decoy-generation") {
        run_decoy_generation(cli_args);
    }
    if let Some(cli_args) = matches.subcommand_matches("spectrum-splitup") {
        run_spectrum_splitup(cli_args)
    }
    if let Some(cli_args) = matches.subcommand_matches("amino-acid-substitution") {
        run_amino_acid_substitution(cli_args)
    }
    if let Some(cli_args) = matches.subcommand_matches("identification") {
        let ident_args = IdentificationArguments::from_cli_args(cli_args);
        identification_task(&ident_args);
    }
}

