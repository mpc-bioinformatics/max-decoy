extern crate postgres;
extern crate csv;

use std::hash::{Hash, Hasher};
use std::path::Path;

use self::postgres::Connection;

use proteomic::models::mass;
use proteomic::models::persistable::Persistable;

#[derive(Debug, Copy, Clone, PartialEq)]
pub enum ModificationPosition {
    Anywhere,
    NTerminus,
    CTerminus
}

pub struct Modification {
    id: i64,                                // BIGSERIAL
    accession: String,                      // TEXT
    name: String,                           // VARCHAR(255)
    position: ModificationPosition,         // SMALLINT
    is_fix: bool,                           // BOOL                 
    amino_acid_one_letter_code: char,       // CHAR(1)
    mono_mass: i64                          // BIGINT
}

impl Modification {
    pub fn new(accession: &str, name: &str, position: ModificationPosition, is_fix: bool, amino_acid_one_letter_code: char, mono_mass: f64) -> Self {
        return Self{
            id: 0,
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
            Self::position_from_char(position_char),
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

    pub fn to_string(&self) -> String {
        return format!(
            "proteomic::models::amino_acids::Modification\n\taccession => {}\n\tname => {}\n\tposition => {}\n\tis fix => {}\n\tamino acid => {}\n\tmono_mass => {}",
            self.accession,
            self.name,
            Self::position_to_string(self.position),
            self.is_fix,
            self.amino_acid_one_letter_code,
            mass::convert_mass_to_float(self.mono_mass)
        );
    }

    // returns the ascii-decimal-code for the positions first char
    fn position_to_int(position: ModificationPosition) -> i16 {
        return match position {
            ModificationPosition::Anywhere => 65,
            ModificationPosition::CTerminus => 67,
            ModificationPosition::NTerminus => 78
        }
    }

    fn position_to_string(position: ModificationPosition) -> String {
        return match position {
            ModificationPosition::Anywhere => "Anywhere".to_string(),
            ModificationPosition::CTerminus => "C-Terminus".to_string(),
            ModificationPosition::NTerminus => "N-Terminus".to_string()
        }
    }

    // argument is positions first char as ascii-decimal-code
    fn position_from_int(position_int: i16) -> ModificationPosition {
        return match position_int {
            65 => ModificationPosition::Anywhere,
            67 => ModificationPosition::CTerminus,
            78 => ModificationPosition::NTerminus,
            _ => panic!("Could not parse int to ModificationPosititon")
        }
    }

    fn position_from_char(position_char: char) -> ModificationPosition {
        return match position_char.to_ascii_uppercase() {
            'A' => ModificationPosition::Anywhere,
            'C' => ModificationPosition::CTerminus,
            'N' => ModificationPosition::NTerminus,
            _ => panic!("Could not parse int to ModificationPosititon")
        }
    }

    pub fn create_from_csv_file(modification_csv_file_path: &str) -> Box<Vec<Modification>> {
        let mut modifications: Vec<Modification> = Vec::new();
        let csv_path = Path::new(modification_csv_file_path);
        let mut reader = match csv::Reader::from_path(&csv_path) {
            Ok(reader) => reader,
            Err(_) => panic!("something went wrong when reading the modifiication csv file. is the file existing and do you have permission to read it?")
        };
        for row in reader.records() {
            let row = row.unwrap();
            modifications.push(Self::new_from_csv_row(&row));
        }
        return Box::new(modifications);
    }
}

impl Clone for Modification {
    fn clone(&self) -> Modification {
        return Self{
            id: self.id,
            accession: self.accession.clone(),
            name: self.name.clone(),
            position: self.position.clone(),
            is_fix: self.is_fix,
            amino_acid_one_letter_code: self.amino_acid_one_letter_code,
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
        & (self.amino_acid_one_letter_code == other.get_amino_acid_one_letter_code())
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

impl Persistable<Modification, i64, String> for Modification {
    fn get_primary_key(&self) -> i64 {
        return self.id;
    }

    fn find(conn: &Connection, primary_key: &i64) -> Result<Self, String> {
        match conn.query("SELECT * FROM amino_acid_modifications WHERE id = $1 LIMIT 1", &[primary_key]) {
            Ok(rows) =>{
                if rows.len() > 0 {
                    Ok(
                        Modification{
                            id: rows.get(0).get(0),
                            accession: rows.get(0).get(1),
                            name: rows.get(0).get(2),
                            position: Self::position_from_int(rows.get(0).get::<usize, i16>(3)),
                            is_fix: rows.get(0).get(4),
                            amino_acid_one_letter_code: (rows.get(0).get::<usize, String>(5)).chars().next().unwrap(),
                            mono_mass: rows.get(0).get(6)
                        }
                    )
                } else {
                    Err("NOHIT".to_owned())
                }
            },
            Err(err) => Err(err.code().unwrap().code().to_owned())
        }
    }

    fn find_by_unique_identifier(conn: &Connection, unique_identifier: &String) -> Result<Self, String> {
        match conn.query(
            "SELECT * FROM amino_acid_modifications WHERE accession = $1 LIMIT 1",
            &[unique_identifier]
        ) {
            Ok(ref rows) if rows.len() > 0 => Ok(
                Modification{
                    id: rows.get(0).get(0),
                    accession: rows.get(0).get(1),
                    name: rows.get(0).get(2),
                    position: Self::position_from_int(rows.get(0).get::<usize, i16>(3)),
                    is_fix: rows.get(0).get(4),
                    amino_acid_one_letter_code: (rows.get(0).get::<usize, String>(5)).chars().next().unwrap(),
                    mono_mass: rows.get(0).get(6)
                }
            ),
            Ok(_rows) => Err("NOHIT".to_owned()),
            Err(err) => Err(err.code().unwrap().code().to_owned())
        }
    }


    fn create(&mut self, conn: &postgres::Connection) -> Result<(), String> {
        match conn.query(
            Self::get_insert_query(),
            &[&self.accession, &self.name, &Self::position_to_int(self.position), &self.is_fix, &self.amino_acid_one_letter_code.to_string(), &self.mono_mass]
        ) {
            Ok(ref rows) if rows.len() > 0 => {
                self.id =  rows.get(0).get(0);
                return Ok(());
            },
            Ok(_rows) => {
                // zero rows means there are a conflict on insert, so the modification exists already
                match Self::find_by_unique_identifier(conn, &self.accession) {
                    Ok(modification) => {
                        self.id = modification.get_primary_key();
                        return Ok(());
                    },
                    Err(err) => Err(format!("cannot insert nor find modification '{}'\n\toriginal error: {}", self.accession, err))
                }
            }
            Err(err) => Err(err.code().unwrap().code().to_owned())
        }
    }

    fn update(&mut self, conn: &postgres::Connection) -> Result<(), String> {
        match conn.query(
            Self::get_update_query(),
            &[&self.accession, &self.name, &Self::position_to_int(self.position), &self.is_fix, &self.amino_acid_one_letter_code.to_string(), &self.mono_mass]
        ) {
            Ok(ref rows) if rows.len() > 0 => Ok(()),
            Ok(_rows) => Err("NORET".to_owned()),
            Err(err) => Err(err.code().unwrap().code().to_owned())
        }
    }

    fn save(&mut self, conn: &postgres::Connection) -> Result<(), String> {
        if self.is_persisted() {
            return self.update(conn);
        } else {
            return self.create(conn);
        }
    }


    fn get_select_primary_key_by_unique_identifier_query() -> &'static str {
        return "SELECT id FROM amino_acid_modifications WHERE accession = $1 LIMIT 1";
    }

    fn get_insert_query() -> &'static str {
        return "INSERT INTO amino_acid_modifications (accession, name, position, is_fix, amino_acid_one_letter_code, mono_mass) VALUES ($1, $2, $3, $4, $5, $6) ON CONFLICT (accession) DO NOTHING RETURNING id";
    }

    fn get_update_query() -> &'static str {
        return "UPDATE amino_acid_modifications SET accession = $2, name = $3, position = $4, is_fix = $5, amino_acid_one_letter_code = $6, mono_mass = $7 WHERE id = $1";
    }

    fn exec_select_primary_key_by_unique_identifier_statement(&mut self, prepared_statement: &postgres::stmt::Statement) -> Result<(), String> {
        match prepared_statement.query(&[&self.accession]) {
            Ok(ref rows) if rows.len() > 0 => {
                self.id = rows.get(0).get(0);
                return Ok(());
            },
            Ok(_rows) => Err("NOHIT".to_owned()),
            Err(err) => Err(err.code().unwrap().code().to_owned())
        }
    }

    fn exec_insert_statement(&mut self, prepared_statement: &postgres::stmt::Statement) -> Result<(), String> {
        match prepared_statement.query(&[&self.accession, &self.name, &Self::position_to_int(self.position), &self.is_fix, &self.amino_acid_one_letter_code.to_string(), &self.mono_mass]) {
            Ok(ref rows) if rows.len() > 0 => {
                self.id = rows.get(0).get(0);
                return Ok(());
            },
            Ok(_rows) => Err("NORET".to_owned()),
            Err(err) => Err(err.code().unwrap().code().to_owned())
        }
    }

    fn exec_update_statement(&mut self, prepared_statement: &postgres::stmt::Statement) -> Result<(), String> {
        match prepared_statement.query(&[&self.accession, &self.name, &Self::position_to_int(self.position), &self.is_fix, &self.amino_acid_one_letter_code.to_string(), &self.mono_mass]) {
            Ok(ref rows) if rows.len() > 0 => Ok(()),
            Ok(_rows) => Err("NORET".to_owned()),
            Err(err) => Err(err.code().unwrap().code().to_owned())
        }
    }


    fn is_persisted(&self) -> bool {
        return self.id < 1;
    }

    fn get_count(conn: &postgres::Connection) -> i64 {
        return match conn.query("SELECT cast(count(id) AS BIGINT) FROM amino_acid_modifications", &[]) {
            Ok(ref rows) if rows.len() > 0 => rows.get(0).get::<usize, i64>(0),
            _ => -1
        };
    }

}