extern crate postgres;
extern crate csv;

use std::hash::{Hash, Hasher};
use std::path::Path;

use proteomic::models::mass;
use proteomic::models::persistable::{handle_postgres_error, Persistable, QueryError, FromSqlRowError};
use proteomic::utility::database_connection::DatabaseConnection;

// it is very important that this dummy is in inserted into the databse first. 
// because postgresql does not consider null-values during testing of unique constraints the follwoing two sql commands to create a modified decoy will both succeded_
//     1. insert into modified_decoys (base_decoy_id, c_terminus_modification_id, n_terminus_modification_id, modification_ids, weight) values (204, null, null, '{null, null, null, null}', 1000000);
//     2. insert into modified_decoys (base_decoy_id, c_terminus_modification_id, n_terminus_modification_id, modification_ids, weight) values (204, null, null, '{null, null, null, null}', 1000000);
// in fact we ended with two redundant ModifiedDecoys, which is bad.
// the plan is to use the id of the dummy modification instead of allowing null during insert.
// vice versa we can replace the dummy modification with None during parsing sql results to object.
// id = -1 will not interfere with any other modification, because BIGSERIAL begins per default with 1.
const DUMMY_MODIFICATION_VALUES: (i64, &str, &str, ModificationPosition, bool, char, i64) = (-1, "max-decoy-internal:dummy", "DO NOT DELETE THIS MODIFICATION!", ModificationPosition::Anywhere, true, '!', 0);


#[derive(Debug, Copy, Clone, PartialEq)]
pub enum ModificationPosition {
    Anywhere,
    NTerminus,
    CTerminus
}

impl ModificationPosition {
    fn to_string(&self) -> String {
        return match self {
            ModificationPosition::Anywhere => "Anywhere".to_string(),
            ModificationPosition::CTerminus => "C-Terminus".to_string(),
            ModificationPosition::NTerminus => "N-Terminus".to_string()
        }
    }

    // returns the ascii-decimal-code for the positions first char
    fn as_small_int(&self) -> i16 {
        return match self {
            ModificationPosition::Anywhere => 65,
            ModificationPosition::CTerminus => 67,
            ModificationPosition::NTerminus => 78
        }
    }

    // argument is positions first char as ascii-decimal-code
    fn from_small_int(position_int: i16) -> ModificationPosition {
        return match position_int {
            65 => ModificationPosition::Anywhere,
            67 => ModificationPosition::CTerminus,
            78 => ModificationPosition::NTerminus,
            _ => panic!("proteomic::models::amino_acids::modification::ModificationPosition::from_small_int(): Could not parse int16 to ModificationPosititon")
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
    id: i64,                                // BIGSERIAL
    accession: String,                      // TEXT
    name: String,                           // VARCHAR(255)
    position_as_int: i16,                   // SMALLINT
    position: ModificationPosition,         
    is_fix: bool,                           // BOOL                 
    amino_acid_one_letter_code: String,     // CHAR(1), must be handled as String, because char does not implemen postgres::types::ToSql, so we can not pass references to .query()
    mono_mass: i64                          // BIGINT
}

impl Modification {
    pub fn new(accession: &str, name: &str, position: ModificationPosition, is_fix: bool, amino_acid_one_letter_code: char, mono_mass: f64) -> Self {
        return Self{
            id: 0,
            accession: accession.to_owned().to_lowercase().trim().to_string(),
            name: name.to_owned().trim().to_string(),
            position_as_int: position.as_small_int(),
            position: position,
            is_fix: is_fix,
            amino_acid_one_letter_code: amino_acid_one_letter_code.to_ascii_uppercase().to_string(),
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
        return match self.amino_acid_one_letter_code.chars().next() {
            Some(one_letter_code) => one_letter_code,
            None => '_'
        };
    }

    pub fn to_string(&self) -> String {
        return format!(
            "proteomic::models::amino_acids::modification::Modification\n\tid => {}\n\taccession => {}\n\tname => {}\n\tposition => {}\n\tis fix => {}\n\tamino acid => {}\n\tmono_mass => {}",
            self.id,
            self.accession,
            self.name,
            self.position.to_string(),
            self.is_fix,
            self.amino_acid_one_letter_code,
            mass::convert_mass_to_float(self.mono_mass)
        );
    }

    pub fn create_from_csv_file(modification_csv_file_path: &str) -> Box<Vec<Modification>> {
        let conn: postgres::Connection = DatabaseConnection::get_database_connection();
        let mut modifications: Vec<Modification> = Vec::new();
        let csv_path = Path::new(modification_csv_file_path);
        let mut reader = match csv::Reader::from_path(&csv_path) {
            Ok(reader) => reader,
            Err(_) => panic!("something went wrong when reading the modifiication csv file. is the file existing and do you have permission to read it?")
        };
        for row in reader.records() {
            let row = row.unwrap();
            modifications.push(Self::new_from_csv_row(&row));
            if let Some(modification) = modifications.last_mut() {
                modification.create(&conn);
            };
        }
        return Box::new(modifications);
    }

    pub fn get_dummy_modification() -> Self {
        return Self {
            id: DUMMY_MODIFICATION_VALUES.0,
            accession: DUMMY_MODIFICATION_VALUES.1.to_owned(),
            name: DUMMY_MODIFICATION_VALUES.2.to_owned(),
            position_as_int: DUMMY_MODIFICATION_VALUES.3.as_small_int(),
            position: DUMMY_MODIFICATION_VALUES.3,
            is_fix: DUMMY_MODIFICATION_VALUES.4,
            amino_acid_one_letter_code: DUMMY_MODIFICATION_VALUES.5.to_string(),
            mono_mass: DUMMY_MODIFICATION_VALUES.6
        };
    }

    pub fn make_sure_dummy_modification_exists() {
        let dummy_modification = Self::get_dummy_modification();
        let conn = DatabaseConnection::get_database_connection();
        let mut query: String = format!("SELECT EXISTS (SELECT id FROM {} WHERE id = $1)", Self::get_table_name());
        let dummmy_modification_exists = match conn.query(query.as_str(), &[&dummy_modification.get_primary_key()]) {
            Ok(ref rows) if rows.len() > 0 => rows.get(0).get::<usize, bool>(0),
            Ok(_rows) => panic!("Panic [proteomic::models::amino_acid::modification::Modification::make_sure_dummy_modification_exists()]: Could not check if dummy-modification exists. Got not return from sql server."),
            Err(err) => panic!("Panic [proteomic::models::amino_acid::modification::Modification::make_sure_dummy_modification_exists()]: {}", handle_postgres_error(&err))
        };
        query = format!("INSERT INTO {} (id, accession, name, position, is_fix, amino_acid_one_letter_code, mono_mass) VALUES ($1, $2, $3, $4, $5, $6, $7)", Self::get_table_name());
        if !dummmy_modification_exists {
            match conn.execute(query.as_str(), &[&dummy_modification.get_primary_key(), &dummy_modification.get_accession(), &dummy_modification.get_name(), &dummy_modification.get_position().as_small_int(), &dummy_modification.is_fix(), &dummy_modification.get_amino_acid_one_letter_code().to_string(), &dummy_modification.get_mono_mass()]) {
                Ok(_) => (),
                Err(err) => panic!("Panic [proteomic::models::amino_acid::modification::Modification::make_sure_dummy_modification_exists()]: {}", handle_postgres_error(&err))
            }
        }
    }
}

impl Clone for Modification {
    fn clone(&self) -> Modification {
        return Self{
            id: self.id,
            accession: self.accession.clone(),
            name: self.name.clone(),
            position_as_int: self.position_as_int,
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

impl Persistable<Modification, i64, String> for Modification {
    fn from_sql_row(row: &postgres::rows::Row) -> Result<Self, FromSqlRowError> {
        let position_as_int: i16 = row.get(3);
        return Ok (
            Self {
                id: row.get(0),
                accession: row.get(1),
                name: row.get(2),
                position_as_int: position_as_int,
                position: ModificationPosition::from_small_int(position_as_int),
                is_fix: row.get(4),
                amino_acid_one_letter_code: row.get::<usize, String>(5),
                mono_mass: row.get(6)
            }
        )
    }

    fn set_primary_key_from_sql_row(&mut self, row: &postgres::rows::Row) {
        self.id = row.get(0);
    }

    fn invalidate_primary_key(&mut self) {
        self.id = 0;
    }

    fn get_primary_key(&self) -> i64 {
        return self.id;
    }

    fn get_table_name() -> &'static str {
        return "amino_acid_modifications";
    }

    fn is_persisted(&self) -> bool {
        return self.id < 1;
    }

    fn find_query() -> &'static str {
        return "SELECT * FROM amino_acid_modifications WHERE id = $1 LIMIT 1;";
    }

    fn create_query() -> &'static str {
        return "INSERT INTO amino_acid_modifications (accession, name, position, is_fix, amino_acid_one_letter_code, mono_mass) VALUES ($1, $2, $3, $4, $5, $6) ON CONFLICT (accession) DO NOTHING RETURNING id;";
    }

    fn create_attributes(&self) -> Box<Vec<&postgres::types::ToSql>>{
        return Box::new(vec![&self.accession, &self.name, &self.position_as_int, &self.is_fix, &self.amino_acid_one_letter_code, &self.mono_mass]);
    }

    fn update_query() -> &'static str{
        return "UPDATE amino_acid_modifications SET accession = $2, name = $3, position = $4, is_fix = $5, amino_acid_one_letter_code = $6, mono_mass = $7 WHERE id = $1;";
    }

    fn update_attributes(&self) -> Box<Vec<&postgres::types::ToSql>>{
        return Box::new(vec![&self.id, &self.accession, &self.name, &self.position_as_int, &self.is_fix, &self.amino_acid_one_letter_code, &self.mono_mass])
    }

    fn delete_query() -> &'static str {
        return "DELETE FROM amino_acid_modifications WHERE id = $1;";
    }

    fn delete_attributes(&self) -> Box<Vec<&postgres::types::ToSql>> {
        return Box::new(vec![&self.id]);
    }

    fn delete_all_query() -> &'static str {
        return "DELETE FROM amino_acid_modifications WHERE id IS NOT NULL;";
    }

    fn exists_query() -> &'static str {
        return "SELECT id FROM amino_acid_modifications WHERE accession = $1 LIMIT 1;";
    }

    fn exists_attributes(&self) -> Box<Vec<&postgres::types::ToSql>> {
        return Box::new(vec![&self.accession]);
    }

    fn before_delete_hook(&self) -> Result<(), QueryError> {return Ok(());}
}