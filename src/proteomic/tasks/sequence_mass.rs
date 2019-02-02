use proteomic::models::mass;
use proteomic::models::amino_acids::amino_acid::AminoAcid;

pub struct SequenceMassArguments {
    sequence: String
}

impl SequenceMassArguments {
    pub fn get_sequence(&self) -> &str {
        return self.sequence.as_str();
    }

    pub fn from_cli_args(cli_args: &clap::ArgMatches) -> Self {
        let sequence: &str = match cli_args.value_of("SEQUENCE") {
            Some(sequence) => sequence,
            None => panic!("proteomic::tasks::digestion::SequenceMassArguments.from_cli_args(): No sequence spezified.")
        };
        return Self {
            sequence: sequence.to_owned()
        }
    }
}

pub fn sequence_mass_task(sequence_mass_args: &SequenceMassArguments) {
    let mass = AminoAcid::get_sequence_weight(sequence_mass_args.get_sequence());
    println!("sequence '{}', has a mass of {} Da", sequence_mass_args.get_sequence(), mass::convert_mass_to_float(mass));
}