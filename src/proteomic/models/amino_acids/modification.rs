use std::hash::{Hash, Hasher};
use std::path::Path;

use proteomic::models::mass;
use proteomic::utility::database_connection::DatabaseConnection;
use proteomic::models::amino_acids::amino_acid::AminoAcid;


#[derive(Debug, Copy, Clone, PartialEq)]
pub enum ModificationPosition {
    Anywhere,
    NTerminus,
    CTerminus
}

impl ModificationPosition {
    #[allow(dead_code)]
    fn to_string(&self) -> String {
        return match self {
            ModificationPosition::Anywhere => "Anywhere".to_string(),
            ModificationPosition::CTerminus => "C-Terminus".to_string(),
            ModificationPosition::NTerminus => "N-Terminus".to_string()
        }
    }

    fn from_char(position_char: char) -> ModificationPosition {
        return match position_char.to_ascii_uppercase() {
            'A' => ModificationPosition::Anywhere,
            'C' => ModificationPosition::CTerminus,
            'N' => ModificationPosition::NTerminus,
            _ => panic!("proteomic::models::amino_acids::modification::ModificationPosition::from_char(): Could not parse char to ModificationPosititon")
        }
    }
}

pub struct Modification {
    accession: String,
    name: String,
    position: ModificationPosition,
    is_fix: bool,
    amino_acid_one_letter_code: char,
    mono_mass: i64
}

impl Modification {
    pub fn new(accession: &str, name: &str, position: ModificationPosition, is_fix: bool, amino_acid_one_letter_code: char, mono_mass: f64) -> Self {
        return Self{
            accession: accession.to_owned().to_lowercase().trim().to_string(),
            name: name.to_owned().trim().to_string(),
            position: position,
            is_fix: is_fix,
            amino_acid_one_letter_code: amino_acid_one_letter_code.to_ascii_uppercase(),
            mono_mass: mass::convert_mass_to_int(mono_mass),
        };
    }

    fn new_from_csv_row(row: &csv::StringRecord) -> Modification {
        // check if row has correct length
        if row.len() != 6 {
            panic!("row has wrong length");
        }
        // type casts
        let position_char: char = row[2].trim().chars().next().expect("can not get position char");
        let is_fix_int: i8 = row[3].trim().parse::<i8>().expect("cannot parse is fix to integer");
        let amino_acid_one_letter_code: char = row[4].trim().chars().next().expect("can not get amino acid one letter code");
        let mono_mass: f64 = row[5].trim().parse::<f64>().expect("can not parse mono mass to float");
        // return modification
        return Self::new(
            &row[0],
            &row[1],
            ModificationPosition::from_char(position_char),
            is_fix_int > 0,
            amino_acid_one_letter_code,
            mono_mass
        )
    }

    pub fn get_mono_mass(&self) -> i64 {
        return self.mono_mass;
    }

    pub fn get_accession(&self) -> &str {
        return &self.accession;
    }

    pub fn get_name(&self) -> &str {
        return &self.name;
    }

    pub fn get_position(&self) -> ModificationPosition {
        return self.position;
    }

    pub fn is_fix(&self) -> bool {
        return self.is_fix;
    }

    pub fn get_amino_acid_one_letter_code(&self) -> char {
        return self.amino_acid_one_letter_code;
    }

    #[allow(dead_code)]
    pub fn to_string(&self) -> String {
        return format!(
            "proteomic::models::amino_acids::modification::Modification\n\taccession => {}\n\tname => {}\n\tposition => {}\n\tis fix => {}\n\tamino acid => {}\n\tmono_mass => {}",
            self.accession,
            self.name,
            self.position.to_string(),
            self.is_fix,
            self.amino_acid_one_letter_code,
            mass::convert_mass_to_float(self.mono_mass)
        );
    }

    pub fn create_from_csv_file(modification_csv_file_path: &str) -> Box<Vec<Modification>> {
        let mut modifications: Vec<Modification> = Vec::new();
        let csv_path = Path::new(modification_csv_file_path);
        let mut reader = match csv::Reader::from_path(&csv_path) {
            Ok(reader) => reader,
            Err(_) => panic!("something went wrong when reading the modifiication csv file. is the file existing and do you have permission to read it?")
        };
        for row_result in reader.records() {
            let row = match row_result {
                Ok(row) => row,
                Err(err) => panic!("proteomic::models::amino_acids::modification::Modification::create_from_csv_file(): Error reading csv-line, see: {}", err)
            };
            modifications.push(Self::new_from_csv_row(&row));
        }
        return Box::new(modifications);
    }

    pub fn to_comet_static_modification_param(&self) -> String {;
        let amino_acid = AminoAcid::get(self.amino_acid_one_letter_code);
        if self.amino_acid_one_letter_code.to_ascii_uppercase() != 'J' {
            return  format!(
                "add_{}_{} = {}",
                self.amino_acid_one_letter_code.to_ascii_uppercase(),
                amino_acid.get_name().to_lowercase(),
                mass::convert_mass_to_float(self.mono_mass)
            );
        } else {
            return format!(
                "add_{}_user_amino_acid = {}",
                self.amino_acid_one_letter_code.to_ascii_uppercase(),
                mass::convert_mass_to_float(amino_acid.get_mono_mass() + self.mono_mass)
            );
        }
    }

    pub fn to_comet_variable_modification_param(&self, modification_number: u8, max_number_of_variable_modification_per_peptide: u8) -> String {
        if modification_number > 9 { panic!("proteomic::models::amino_acids::modification::Modification.to_comet_variable_modification_param(): modification_number is not a number from 0 to 9") }
        let distant_to_terminus: i8 = match self.position {
            ModificationPosition::Anywhere => -1,
            ModificationPosition::CTerminus => 1,
            ModificationPosition::NTerminus => 1
        };
        let distant_refere_to_terminus: u8 = match self.position {
            ModificationPosition::Anywhere => 2,
            ModificationPosition::CTerminus => 3,
            ModificationPosition::NTerminus => 2
        };
        return format!(
            "variable_mod0{} = {} {} 0 {} {} {} 0",
            modification_number,
            mass::convert_mass_to_float(self.mono_mass),
            self.get_amino_acid_one_letter_code().to_uppercase(),
            max_number_of_variable_modification_per_peptide,
            distant_to_terminus,
            distant_refere_to_terminus
        );
    }
}

impl Clone for Modification {
    fn clone(&self) -> Modification {
        return Self{
            accession: self.accession.clone(),
            name: self.name.clone(),
            position: self.position.clone(),
            is_fix: self.is_fix,
            amino_acid_one_letter_code: self.amino_acid_one_letter_code.clone(),
            mono_mass: self.mono_mass
        };
    }
}

// PartialEq-implementation to use this type in a HashSet
// modifcation name is not considered here because the modification might have an additional name or spelling
impl PartialEq for Modification {
    fn eq(&self, other: &Modification) -> bool {
       return (self.accession == other.get_accession())
        & (self.position == other.get_position())
        & (self.is_fix == other.is_fix())
        & (self.get_amino_acid_one_letter_code() == other.get_amino_acid_one_letter_code())
        & (self.mono_mass == other.get_mono_mass());
    }
}

// Eq-implementation to use this type in a HashSet
impl Eq for Modification {}

// Hash-implementation to use this type in a HashSet
impl Hash for Modification {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.accession.hash(state);
    }
}