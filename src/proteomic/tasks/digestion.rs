use proteomic::utility::input_file_digester::file_digester::FileDigester;
use proteomic::utility::input_file_digester::fasta_digester::FastaDigester;

pub struct DigestionArguments {
    input_file: String,
    input_format: String,
    thread_count: usize,
    number_of_missed_cleavages: u8,
    min_peptide_length: usize,
    max_peptide_length: usize,
    enzym_name: String
}

impl DigestionArguments {
    pub fn get_input_file(&self) -> &str {
        return self.input_file.as_str();
    }

    pub fn get_input_format(&self) -> &str {
        return self.input_format.as_str();
    }

    pub fn get_thread_count(&self) -> usize {
        return self.thread_count;
    }

    pub fn get_number_of_missed_cleavages(&self) -> u8 {
        return self.number_of_missed_cleavages;
    }

    pub fn get_min_peptide_length(&self) -> usize {
        return self.min_peptide_length;
    }

    pub fn get_max_peptide_length(&self) -> usize {
        return self.max_peptide_length;
    }

    pub fn get_enzym_name(&self) -> &str {
        return self.enzym_name.as_str();
    }

    pub fn from_cli_args(cli_args: &clap::ArgMatches) -> Self {
        let input_file: &str = match cli_args.value_of("INPUT_FILE") {
            Some(file) => file,
            None => panic!("proteomic::tasks::digestion::DigestionToFileArguments.from_cli_args(): No input-file spezified.")
        };
        let input_format: &str = match cli_args.value_of("INPUT_FORMAT") {
            Some(format) => format,
            None => panic!("proteomic::tasks::digestion::DigestionToFileArguments.from_cli_args(): No input-format spezified.")
        };
        let thread_count: usize = match cli_args.value_of("THREAD_COUNT") {
            Some(number_string) => match number_string.to_owned().parse::<usize>() {
                Ok(mut count) => count,
                Err(_err) => panic!("proteomic::tasks::digestion::DigestionToFileArguments.from_cli_args(): Could not parse thread-count to unsigned integer."),
            },
            None => {
                println!("WARNING [proteomic::tasks::digestion::DigestionToFileArguments.from_cli_args()]: No thread-count spezified, set it to 1.");
                1
            }
        };
        let number_of_missed_cleavages: u8 = match cli_args.value_of("NUMBER_OF_MISSED_CLEAVAGES") {
            Some(count) => {
                match count.to_owned().parse::<u8>() {
                    Ok(parsed_count) => parsed_count,
                    Err(_err) => panic!("proteomic::tasks::digestion::DigestionToFileArguments.from_cli_args(): Could not parse number-of-missed-cleavages to unsigned integer.")
                }
            },
            None => {
                println!("WARNING [proteomic::tasks::digestion::DigestionToFileArguments.from_cli_args()]: number-of-missed-cleavages not spezified, set it to 0.");
                0
            }
        };
        if number_of_missed_cleavages > 60 {
            panic!("proteomic::tasks::digestion::DigestionToFileArguments.from_cli_args(): Possible maximal value for number-of-missed-cleavages is 60.");
        }
        let min_peptide_length: usize = match cli_args.value_of("MIN_PEPTIDE_LENGTH") {
            Some(count) => {
                match count.to_owned().parse::<usize>() {
                    Ok(parsed_count) => parsed_count,
                    Err(_err) => panic!("proteomic::tasks::digestion::DigestionToFileArguments.from_cli_args(): Could not parse min-peptide-length to unsigned integer.")
                }
            },
            None => {
                println!("WARNING [proteomic::tasks::digestion::DigestionToFileArguments.from_cli_args()]: min-peptide-length length not spezified, set it to 6");
                6
            }
        };
        let max_peptide_length: usize = match cli_args.value_of("MAX_PEPTIDE_LENGTH") {
            Some(count) => {
                match count.to_owned().parse::<usize>() {
                    Ok(parsed_count) => parsed_count,
                    Err(_err) =>  panic!("proteomic::tasks::digestion::DigestionToFileArguments.from_cli_args(): Could not parse max-peptide-length to unsigned integer.")
                }
            },
            None => {
                println!("WARNING [proteomic::tasks::digestion::DigestionToFileArguments.from_cli_args()]: max-peptide-length length not spezified, set it to 50");
                50
            }
        };
        if max_peptide_length > 60 {
            panic!("proteomic::tasks::digestion::DigestionToFileArguments.from_cli_args(): Possible maximal value for max-peptide-length is 60.");
        }
        let enzym_name: &str = match cli_args.value_of("ENZYM_NAME") {
            Some(enzym) => enzym,
            None => {
                println!("WARNING [proteomic::tasks::digestion::DigestionToFileArguments.from_cli_args()]: No enzym spezified, use Trypsin.");
                "trypsin"
            }
        };
        if min_peptide_length > max_peptide_length {
            panic!("proteomic::tasks::digestion::DigestionToFileArguments.from_cli_args(): min-peptide-length must be less or equals than max-peptide-length");
        }
        return Self {
            input_file: input_file.to_owned(),
            input_format: input_format.to_owned(),
            thread_count: thread_count,
            number_of_missed_cleavages: number_of_missed_cleavages,
            min_peptide_length: min_peptide_length,
            max_peptide_length: max_peptide_length,
            enzym_name: enzym_name.to_owned()
        }
    }
}

pub fn digest_to_database_task(digestion_arguments: &DigestionArguments) {
    let mut digester = FastaDigester::new(
        digestion_arguments.get_input_file(),
        digestion_arguments.get_thread_count(),
        digestion_arguments.get_number_of_missed_cleavages(),
        digestion_arguments.get_min_peptide_length(),
        digestion_arguments.get_max_peptide_length()
    );
    let seconds = digester.process_file(digestion_arguments.get_enzym_name());
    println!("need {} days", seconds / 60.0 / 60.0 / 24.0)
}