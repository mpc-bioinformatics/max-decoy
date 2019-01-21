extern crate clap;
extern crate num_cpus;
extern crate postgres;
extern crate postgres_array;
extern crate quick_xml;
extern crate sha1;

use std::collections::{HashSet, HashMap};

use clap::{Arg, App, SubCommand};

mod proteomic;
use proteomic::utility::input_file_digester::file_digester::FileDigester;
use proteomic::utility::input_file_digester::fasta_digester::FastaDigester;
use proteomic::utility::database_connection::DatabaseConnection;
use proteomic::utility::decoy_generator::DecoyGenerator;
use proteomic::utility::mz_ml::mz_ml_reader::MzMlReader;
use proteomic::utility::mz_ml::spectrum::Spectrum;
use proteomic::utility::mz_ml::chromatogram::Chromatrogram;

use proteomic::models::enzyms::trypsin::Trypsin;
use proteomic::models::persistable::Persistable;
use proteomic::models::peptide::Peptide;
use proteomic::models::decoys::decoy::PlainDecoy;
use proteomic::models::amino_acids::modification::Modification;


fn run_digestion(digest_cli_args: &clap::ArgMatches) {
    let mut error_in_digest_args = false;
    let input_file: &str = match digest_cli_args.value_of("INPUT_FILE") {
        Some(file) => file,
        None => {
            println!("ERROR [digest]: No input file spezified.");
            error_in_digest_args = true;
            ""
        }
    };
    let input_format: &str = match digest_cli_args.value_of("INPUT_FORMAT") {
        Some(format) => format,
        None => {
            println!("ERROR [digest]: No input format spezified.");
            ""
        }
    };
    let cpu_thread_count: usize = num_cpus::get();
    let thread_count: usize = match digest_cli_args.value_of("THREAD_COUNT") {
        Some(count) => {
            match count.to_lowercase().as_str() {
                "max" => cpu_thread_count,
                _ => match count.to_owned().parse::<usize>() {
                    Ok(mut parsed_count) => {
                        let cpu_thread_count: usize = num_cpus::get();
                        if (1 > parsed_count) | (parsed_count > num_cpus::get()) {
                            println!("ERROR [digest]: Threadcount must be between {} and {}.", 1, cpu_thread_count);
                            error_in_digest_args = true;
                        }
                        parsed_count
                    } ,
                    Err(_err) => {
                        println!("ERROR [digest]: Threadcount is not a positive integer.");
                        error_in_digest_args = true;
                        0
                    }
                }
            }
        },
        None => {
            println!("WARNING [digest]: No thread spezified, set it to 1.");
            1
        }
    };
    let number_of_missed_cleavages: i16 = match digest_cli_args.value_of("NUMBER_OF_MISSED_CLEAVAGES") {
        Some(count) => {
            match count.to_owned().parse::<i16>() {
                Ok(parsed_count) => parsed_count,
                Err(_err) => {
                    println!("ERROR [digest]: Number of missed cleavages is not a positive integer.");
                    error_in_digest_args = true;
                    1
                }
            }
        },
        None => {
            println!("WARNING [digest]: Number of missed cleavages not spezified, set it to 2.");
            2
        }
    };
    let min_peptide_length: usize = match digest_cli_args.value_of("MIN_PEPTIDE_LENGTH") {
        Some(count) => {
            match count.to_owned().parse::<usize>() {
                Ok(parsed_count) => parsed_count,
                Err(_err) => {
                    println!("ERROR [digest]: Minimum peptide length is not an positive integer.");
                    error_in_digest_args = true;
                    1
                }
            }
        },
        None => {
            println!("WARNING [digest]: Minimum peptide length not spezified, set it to 6");
            6
        }
    };
    let max_peptide_length: usize = match digest_cli_args.value_of("MAX_PEPTIDE_LENGTH") {
        Some(count) => {
            match count.to_owned().parse::<usize>() {
                Ok(parsed_count) => parsed_count,
                Err(_err) => {
                    println!("ERROR [digest]: Max peptide length is not an positive integer.");
                    error_in_digest_args = true;
                    1
                }
            }
        },
        None => {
            println!("WARNING [digest]: Max peptide length not spezified, set it to 50");
            50
        }
    };
    let enzym_name: &str = match digest_cli_args.value_of("ENZYM_NAME") {
        Some(enzym) => enzym,
        None => {
            println!("WARNING [digest]: No enzym spezified, use Trypsin.");
            "trypsin"
        }
    };
    if !error_in_digest_args {
        let database_connection: postgres::Connection = DatabaseConnection::get_database_connection();
        let mut results_for_digest_and_commit: (usize, usize, f64) = (0, 0, 0.0);
        let mut results_for_counting: (usize, usize) = (0, 0);
        match input_format {
            "fasta" => {
                let mut file_handler = match enzym_name.to_lowercase().as_str() {
                    "trypsin" => FastaDigester::<Trypsin>::new(
                        input_file,
                        thread_count,
                        number_of_missed_cleavages,
                        min_peptide_length,
                        max_peptide_length
                    ),
                    _ => FastaDigester::<Trypsin>::new(
                        input_file,
                        thread_count,
                        number_of_missed_cleavages,
                        min_peptide_length,
                        max_peptide_length
                    )
                };
                results_for_digest_and_commit = file_handler.process_file();
                println!(
                    "{:<20}{:<20}",
                    "type",
                    "comitted"
                );
                println!(
                    "{:<20}{:<20}",
                    "proteins",
                    results_for_digest_and_commit.0
                );
                println!(
                    "{:<20}{:<20}",
                    "peptides",
                    results_for_digest_and_commit.1
                );
                println!(
                    "{:<20}{:<20}",
                    "commit time",
                    results_for_digest_and_commit.2
                );
            },
            _ => println!("ERROR [digest]: Input format unknown.")
        }
    } else {
        println!("ERROR [digest]: Did nothing, some error occured.")
    }
}

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
    let number_of_decoys_per_target: i64 = match decoy_generation_cli_args.value_of("NUMBER_OF_DECOYS") {
        Some(number_string) => {
            match number_string.to_owned().parse::<i64>() {
                Ok(number) => number,
                Err(_) => panic!("ERROR [decoy-generation]: could not cast number-of-decoys to integer")
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
    let mz_ml_file: &str = match decoy_generation_cli_args.value_of("MZ_ML_FILE") {
        Some(mz_ml_file) => mz_ml_file,
        None => panic!("you must provide a MzMl-file")
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
    let mz_ml_reader = MzMlReader::new(mz_ml_file);
    let spectra: Vec<Spectrum> = *mz_ml_reader.get_ms_two_spectra();
    let conn = DatabaseConnection::get_database_connection();
    for spectrum in spectra {
        let generator: DecoyGenerator = DecoyGenerator::new(*spectrum.get_precurso_mass(), upper_mass_tolerance, lower_mass_tolerance, thread_count, max_modifications_per_decoy, &fixed_modifications_map, &variable_modifications_map);
        let targets = match Peptide::find_where(&conn, "weight BETWEEN $1 AND $2", &[&generator.get_lower_weight_limit(), &generator.get_upper_weight_limit()]) {
            Ok(targets) => targets,
            Err(err) => panic!("main::run_decoy_generation() could not gether target count: {}", err)
        };
        println!("found {} targets", targets.len());
        let decoy_count = match PlainDecoy::count_where_mass_tolerance(&conn, generator.get_lower_weight_limit(), generator.get_upper_weight_limit()) {
            Ok(count) => count,
            Err(err) => panic!("main::run_decoy_generation() could not gether decoy count: {}", err)
        };
        let mut number_of_decoys_to_generate: i64 = (targets.len() as i64 * number_of_decoys_per_target) - decoy_count;
        if number_of_decoys_to_generate < 0 {
            number_of_decoys_to_generate = 0;
        }
        let number_of_decoys_to_generate_freeze = number_of_decoys_to_generate;
        println!("need to generate {} decoys", number_of_decoys_to_generate);
        let start_time: f64 = time::precise_time_s();
        if number_of_decoys_to_generate > 0 {
            number_of_decoys_to_generate -= generator.vary_targets(&targets, number_of_decoys_to_generate as usize) as i64;
            println!("still {} decoys needed after shuffling and varying targets", number_of_decoys_to_generate);
        }
        if number_of_decoys_to_generate > 0 {
            number_of_decoys_to_generate -= generator.generate_decoys(number_of_decoys_to_generate as usize) as i64;
        }
        let stop_time: f64 = time::precise_time_s();
        println!("generate {} decoys in {} s", number_of_decoys_to_generate_freeze, stop_time - start_time);
        // get decoys
        // write to fasta files
        // call comet
        break;
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
    let chromatograms: Vec<Chromatrogram> = *mz_ml_reader.get_chromatograms();
    let content_before_spectrum_list = mz_ml_reader.get_content_before_spectrum_list();
    let spectra: Vec<Spectrum> = *mz_ml_reader.get_ms_two_spectra();
    for spectrum in spectra.iter() {
        spectrum.to_mz_ml(content_before_spectrum_list.as_str(), &chromatograms, destination_folder, file_suffix);
    }
    let stop_time: f64 = time::precise_time_s();
    println!("found and write {} MS2-spectra in {} s", spectra.len(), stop_time - start_time);
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
            .value_name("THREADCOUNT")
            .takes_value(true)
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
        )
        .arg(
            Arg::with_name("ENZYM_NAME")
            .short("e")
            .long("enzym-name")
            .value_name("ENZYM_NAME")
            .takes_value(true)
        )
        .arg(
            Arg::with_name("CHECK_DB_VALUES")
            .long("check-db-values")
            .help("Run digestion again, stores peptides in unique list and compares number of peptides in databases with number of peptides in unique list. (ATTENTIONS: Might runs out of RAM, because unique list is stored in RAM)")
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
            Arg::with_name("MZ_ML_FILE")
            .short("z")
            .long("mz-ml-file")
            .value_name("MZ_ML_FILE")
            .takes_value(true)
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
    .get_matches();


    if let Some(cli_args) = matches.subcommand_matches("digest") {
        run_digestion(cli_args);
    }
    if let Some(cli_args) = matches.subcommand_matches("decoy-generation") {
        run_decoy_generation(cli_args);
    }
    if let Some(cli_args) = matches.subcommand_matches("spectrum-splitup") {
        run_spectrum_splitup(cli_args)
    }
}

