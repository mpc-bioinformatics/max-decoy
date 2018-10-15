extern crate clap;
extern crate num_cpus;

use clap::{Arg, App, SubCommand};

mod proteomic;
use proteomic::utility::input_file_digester::file_digester::FileDigester;
use proteomic::utility::input_file_digester::fasta_digester::FastaDigester;
use proteomic::utility::enzyms::trypsin::Trypsin;

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
        match input_format {
            "" => {
                let file_handler = match enzym_name.to_lowercase().as_str() {
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
                file_handler.process_file();
            },
            _ => println!("ERROR [digest]: Input format unknown.")
        }
    } else {
        println!("ERROR [digest]: Did nothing, some error occured.")
    }

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
    )
    .get_matches();


    if let Some(matches) = matches.subcommand_matches("digest") {
        run_digestion(matches);
    }
}