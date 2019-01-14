use std::collections::HashMap;
use std::fmt;

use proteomic::models::amino_acids::modification::{Modification, ModificationPosition};
use proteomic::models::amino_acids::amino_acid::AminoAcid;
use proteomic::models::mass;
use proteomic::models::mass::neutral_loss::NeutralLoss;
use proteomic::models::decoys::decoy::Decoy;
use proteomic::models::decoys::base_decoy::BaseDecoy;
use proteomic::models::decoys::modified_decoy::ModifiedDecoy;
use proteomic::utility::combinations::n_choose_k::NChooseK;
use proteomic::models::persistable::{Persistable, QueryOk, QueryError};
use proteomic::models::peptide::Peptide;

pub enum NewDecoyError {
    ModificationDoesNotMatchToAminoAcid,       // modification does not relate to amino acid
    OutrangeMassTolerance,
    ModificationIsNotFix,
    ModificationIsNotVariable,
    AlreadyFixModificationInPlace
}

impl NewDecoyError {
    pub fn to_string(&self) -> String {
        match self {
            NewDecoyError::ModificationDoesNotMatchToAminoAcid => format!("NewDecoyError::ModificationDoesNotMatchToAminoAcid"),
            NewDecoyError::OutrangeMassTolerance => format!("NewDecoyError::OutrangeMassTolerance"),
            NewDecoyError::ModificationIsNotFix => format!("NewDecoyError::ModificationIsNotFix"),
            NewDecoyError::ModificationIsNotVariable => format!("NewDecoyError::ModificationIsNotVariable"),
            NewDecoyError::AlreadyFixModificationInPlace => format!("NewDecoyError::AlreadyFixModificationInPlace")
        }
    }
}

impl fmt::Display for NewDecoyError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        return write!(f, "{}", self.to_string());
    }
}



pub struct NewDecoy {
    aa_sequence: Vec<char>,
    modifications: Vec<Option<Modification>>,
    weight: i64,
    lower_weight_limit: i64,
    upper_weight_limit: i64,
    n_terminus_modification: Option<Modification>,
    c_terminus_modification: Option<Modification>,
    number_of_modifications: i32
}

impl NewDecoy {
    pub fn new(lower_weight_limit: i64, upper_weight_limit: i64) -> Self {
        let water_loss: NeutralLoss = NeutralLoss::get("H2O");
        return Self {
            aa_sequence: Vec::new(),
            modifications: Vec::new(),
            weight: water_loss.get_mono_mass(),
            lower_weight_limit: lower_weight_limit,
            upper_weight_limit: upper_weight_limit,
            n_terminus_modification: None,
            c_terminus_modification: None,
            number_of_modifications: 0
        }
    }

    pub fn from_string(aa_sequence: &str, lower_weight_limit: i64, upper_weight_limit: i64, fix_modifications: &HashMap<char, Modification>) -> Self {
        let mut new_decoy = Self::new(lower_weight_limit, upper_weight_limit);
        for amino_acids_one_letter_code in aa_sequence.chars() {
            let amino_acid = AminoAcid::get(amino_acids_one_letter_code);
            let modification_option = match fix_modifications.get(&amino_acids_one_letter_code) {
                Some(ref modification) => Some((*modification).clone()), 
                None => None
            };
            match new_decoy.push_amino_acid_and_fix_modification(&amino_acid, &modification_option) {
                Ok(_) => (),
                Err(err) => match err {
                    NewDecoyError::OutrangeMassTolerance => continue,
                    NewDecoyError::ModificationIsNotFix => panic!("proteomic::models::decoys::decoy::new_decoy::NewDecoy::from_string(): variable modification is passed to push_amino_acid_and_fix_modification(), which actually should not happen here."),
                    _ => panic!("proteomic::models::decoys::decoy::new_decoy::NewDecoy::from_string(): only NewDecoyError::OutrangeMassTolerance and NewDecoyError::ModificationIsNotFix should returned here, but something else was retuned."),
                }
            }
        }
        return new_decoy;
    }

    pub fn get_unmodified_weight(&self) -> i64 {
        return AminoAcid::get_sequence_weight(self.get_aa_sequence().as_str());
    }

    pub fn get_number_of_modifications(&self) -> i32 {
        return self.number_of_modifications;
    }

    pub fn hits_mass_tolerance(&self) -> bool {
        return (self.lower_weight_limit <= self.weight) & (self.weight <= self.upper_weight_limit);
    }

    /// returns positive value if weight is less than lower_weight_limit
    /// returns negative value if weight is higher than upper_weight_limit
    /// returns zero if weight is between or equals lower_weight_limit and upper_weight_limit
    pub fn get_distance_to_mass_tolerance(&self) -> i64 {
        if self.weight < self.lower_weight_limit {
            return self.lower_weight_limit - self.weight;
        } else if self.weight > self.upper_weight_limit {
            return self.upper_weight_limit - self.weight;
        } else {
            return 0;
        }
    }

    fn push_amino_acid(&mut self, amino_acid: &AminoAcid) {
        self.weight += amino_acid.get_mono_mass();
        self.aa_sequence.push(amino_acid.get_one_letter_code());
    }

    fn undo_push_amino_acid(&mut self, amino_acid: &AminoAcid) {
        self.weight -= amino_acid.get_mono_mass();
        let last_idx = self.aa_sequence.len() - 1;
        self.aa_sequence.remove(last_idx);
    }

    /// returns the old c_terminus_modification
    fn push_modification(&mut self, modification_option: &Option<Modification>) -> Option<Modification> {
        let old_c_terminus_modification = self.remove_c_terminus_modification();        // after pushing a new amino acid, the old-last one has no c-terminus any more. so remove the c-terminus modification
        if let Some(modification) = modification_option {
            match modification.get_position() {
                ModificationPosition::Anywhere => {
                    self.weight += modification.get_mono_mass();
                    self.number_of_modifications += 1;
                    self.modifications.push(Some(modification.clone()));
                },
                ModificationPosition::CTerminus => {
                    self.weight += modification.get_mono_mass();
                    self.number_of_modifications += 1;
                    self.c_terminus_modification = Some(modification.clone());
                },
                // only if this is the first amino acid and modification
                ModificationPosition::NTerminus if self.aa_sequence.len() == 1 => {
                    self.weight += modification.get_mono_mass();
                    self.number_of_modifications += 1;
                    self.n_terminus_modification = Some(modification.clone());
                },
                _ => self.modifications.push(None)
            }
        } else {
            self.modifications.push(None)
        }
        return old_c_terminus_modification;         // return old c-terminus modification in case of undoing this push
    }

    fn undo_push_modification(&mut self, old_c_terminus_modification_option: &Option<Modification>) {
        // remove current c-terminus-modification
        self.remove_c_terminus_modification();
        // if this is the first amino acid remove the n-terminus modification too
        if self.aa_sequence.len() == 1 {
            self.remove_n_terminus_modification();
        }
        // remove last modification/none of modifications
        let last_idx = self.modifications.len() - 1;
        match self.modifications.get(last_idx) {
            Some(modification_option) => match modification_option {
                Some(modification) => {
                    self.weight -= modification.get_mono_mass();
                    self.number_of_modifications -= 1;
                },
                None => ()
            },
            None => ()
        }
        self.modifications.remove(last_idx);
        self.number_of_modifications -= 1;
        // add old c-terminus modification if the sequence is not empty
        if self.aa_sequence.len() > 0 {
            match old_c_terminus_modification_option {
                Some(modification) => {
                    self.weight += modification.get_mono_mass();
                    self.number_of_modifications += 1;
                    self.c_terminus_modification = Some(modification.clone());
                },
                None => ()
            }
        }
    }

    pub fn as_base_decoy(&self) -> BaseDecoy {
        return BaseDecoy::new(
            self.get_header().as_str(),
            self.get_aa_sequence().as_str(),
            AminoAcid::get_sequence_weight(self.get_aa_sequence().as_str())
        )
    }

    pub fn as_modified_decoy(&self, base_decoy: &BaseDecoy) -> ModifiedDecoy {
        return ModifiedDecoy::new(
            base_decoy,
            self.c_terminus_modification.clone(),
            self.n_terminus_modification.clone(),
            &self.modifications,
            self.weight
        );
    }

    fn remove_c_terminus_modification(&mut self) -> Option<Modification> {
        match self.c_terminus_modification {
            Some(ref modification) => {
                self.weight -= modification.get_mono_mass();
                self.number_of_modifications -= 1;
            },
            None => ()
        }
        let old_c_terminus_modification = self.c_terminus_modification.clone();
        self.c_terminus_modification = None;
        return old_c_terminus_modification;
    }

    fn remove_n_terminus_modification(&mut self) {
        match self.n_terminus_modification {
            Some(ref modification) => {
                self.weight -= modification.get_mono_mass();
                self.number_of_modifications -= 1;
            },
            None => ()
        }
        self.n_terminus_modification = None;
    }

    pub fn push_amino_acid_and_fix_modification(&mut self, amino_acid: &AminoAcid, modification_option: &Option<Modification>) -> Result<(), NewDecoyError> {
        match modification_option {
            Some(modification) => {
                if !modification.is_fix() {
                    return Err(NewDecoyError::ModificationIsNotFix);
                }
                if modification.get_amino_acid_one_letter_code() != amino_acid.get_one_letter_code() {
                    return Err(NewDecoyError::ModificationDoesNotMatchToAminoAcid);
                }
            },
            None => ()
        }
        self.push_amino_acid(amino_acid);
        let old_c_terminus_modifcation = self.push_modification(modification_option);
        if self.get_distance_to_mass_tolerance() < 0 {
            self.undo_push_modification(&old_c_terminus_modifcation);
            self.undo_push_amino_acid(&amino_acid);
            return Err(NewDecoyError::OutrangeMassTolerance);
        }
        return Ok(());
    }

    fn set_variable_c_terminus_modification(&mut self, modification: &Modification) -> Result<(), NewDecoyError> {
        if modification.is_fix() {
            return Err(NewDecoyError::ModificationIsNotVariable);
        }
        if modification.get_amino_acid_one_letter_code() != self.get_c_terminus_amino_acid() {
            return Err(NewDecoyError::ModificationDoesNotMatchToAminoAcid);
        }
        if self.c_terminus_modification.is_some() {
            return Err(NewDecoyError::AlreadyFixModificationInPlace);
        }
        self.weight += modification.get_mono_mass();
        self.number_of_modifications += 1;
        self.c_terminus_modification = Some(modification.clone());
        if self.get_distance_to_mass_tolerance() < 0 {
            return Err(NewDecoyError::OutrangeMassTolerance);
        }
        return Ok(());
    }

    fn set_variable_n_terminus_modification(&mut self, modification: &Modification) -> Result<(), NewDecoyError> {
        if modification.is_fix() {
            return Err(NewDecoyError::ModificationIsNotVariable);
        }
        if modification.get_amino_acid_one_letter_code() != self.get_n_terminus_amino_acid() {
            return Err(NewDecoyError::ModificationDoesNotMatchToAminoAcid);
        }
        if self.n_terminus_modification.is_some() {
            return Err(NewDecoyError::AlreadyFixModificationInPlace);
        }
        self.weight += modification.get_mono_mass();
        self.number_of_modifications += 1;
        self.n_terminus_modification = Some(modification.clone());
        if self.get_distance_to_mass_tolerance() < 0 {
            return Err(NewDecoyError::OutrangeMassTolerance);
        }
        return Ok(());
    }

    pub fn set_variable_modification_at(&mut self, idx: usize, modification: &Modification) -> Result<(), NewDecoyError> {
        if modification.is_fix() {
            return Err(NewDecoyError::ModificationIsNotVariable);
        }      
        if modification.get_amino_acid_one_letter_code() != self.get_amino_acid_at(idx) {
            return Err(NewDecoyError::ModificationDoesNotMatchToAminoAcid);
        }
        if (idx == 0) & (modification.get_position() == ModificationPosition::NTerminus) {
            return self.set_variable_n_terminus_modification(modification);
        } else if (idx == self.aa_sequence.len() - 1) & (modification.get_position() == ModificationPosition::CTerminus) {
            return self.set_variable_c_terminus_modification(modification);
        } else if modification.get_position() == ModificationPosition::Anywhere {
            match self.modifications.get_mut(idx) {
                Some(modification_option) => match modification_option {
                    Some(_) => return Err(NewDecoyError::AlreadyFixModificationInPlace),
                    None => *modification_option = Some(modification.clone())
                },
                None => ()
            }
            self.weight += modification.get_mono_mass();
            self.number_of_modifications += 1;
            if self.get_distance_to_mass_tolerance() < 0 {
                return Err(NewDecoyError::OutrangeMassTolerance);
            }
        }
        return Ok(());
    }

    pub fn remove_all_variable_modifications(&mut self) {
        let mut set_current_modification_to_none = match self.n_terminus_modification {
            Some(ref mut modification) if !modification.is_fix() => {
                self.weight -= modification.get_mono_mass();
                self.number_of_modifications -= 1;
                true
            },
            _ => false
        };
        if set_current_modification_to_none { self.n_terminus_modification = None; }

        set_current_modification_to_none = match self.c_terminus_modification {
            Some(ref mut modification) if !modification.is_fix() => {
                self.weight -= modification.get_mono_mass();
                self.number_of_modifications -= 1;
                true
            },
            _ => false
        };
        if set_current_modification_to_none { self.c_terminus_modification = None; }

        for modification_option in self.modifications.iter_mut() {
            set_current_modification_to_none = match modification_option {
                Some(modification) if !modification.is_fix() => {
                    self.weight -= modification.get_mono_mass();
                    self.number_of_modifications -= 1;
                    true
                },
                _ => false
            };
           if set_current_modification_to_none { *modification_option = None; }
        }
    }

    fn remove_modification_at(&mut self, idx: usize) {
        if idx == 0 { self.remove_n_terminus_modification(); }
        if idx == self.aa_sequence.len() - 1 { self.remove_c_terminus_modification(); }
        match self.modifications.get_mut(idx) {
            Some(modification_options) => {
                match modification_options {
                    Some(ref modification) => {
                        self.weight -= modification.get_mono_mass();
                        self.number_of_modifications -= 1;
                    },
                    None => ()
                }
                *modification_options = None;
            },
            None => ()
        };
    }

    fn add_modification_at(&mut self, idx: usize, modification: &Modification) -> Result<(), NewDecoyError> {
        match self.aa_sequence.get(idx) {
            Some(one_letter_code) => {
                if *one_letter_code != modification.get_amino_acid_one_letter_code() {
                    return Err(NewDecoyError::ModificationDoesNotMatchToAminoAcid);
                }
            },
            None => return Err(NewDecoyError::ModificationDoesNotMatchToAminoAcid)
        }
        if (idx == 0) & (modification.get_position() == ModificationPosition::NTerminus) {
            self.weight += modification.get_mono_mass();
            self.number_of_modifications += 1;
            self.n_terminus_modification = Some(modification.clone());
        } else if (idx == self.aa_sequence.len() - 1) & (modification.get_position() == ModificationPosition::CTerminus) {
            self.weight += modification.get_mono_mass();
            self.number_of_modifications += 1;
            self.c_terminus_modification = Some(modification.clone());
        } else if modification.get_position() == ModificationPosition::Anywhere {
            self.weight += modification.get_mono_mass();
            self.number_of_modifications += 1;
            match self.modifications.get_mut(idx) {
                Some(current_modification) => *current_modification = Some(modification.clone()),
                None => panic!("proteomic::models::decoys::new_decoy::NewDecoy.add_modification_at(): expect Some(Modification) instead of None")
            }
        }
        return Ok(());
    }

    pub fn swap_amino_acids_to_hit_mass_tolerance(&mut self, conn: &postgres::Connection, interchange_map: &HashMap<char, HashMap<char, i64>>, fix_modifications_map: &HashMap<char, Modification>, max_number_of_modifications: u8, varibale_modification_map: &HashMap<char, Modification>) -> usize {
        let mut created_counter: usize = 0;
        let mut tries: u16 = 0;
        'tries: loop {
            'sequence: for idx in 0..self.aa_sequence.len() {
                let increase_weight: bool = match self.get_distance_to_mass_tolerance() {
                    dist if dist > 0 => true,
                    dist if dist < 0 => false,
                    _ => return 1                 // end swap_amino_acids_to_hit_mass_tolerance() here because self already hits mass tolerance
                };
                let aa_one_letter_code = self.get_amino_acid_at(idx);
                if let Some(ref replacements) = interchange_map.get(&aa_one_letter_code) {
                    'replacements: for (aa_replacement, weight_change) in replacements.iter() {
                        // continue if weight must increase but weight_change is negative or 
                        // weight mus decrease but weight_change is positive
                        if increase_weight & (*weight_change < 0) | !increase_weight & (*weight_change > 0) { continue 'replacements; }
                        self.remove_modification_at(idx);
                        // do not use weight_change, because it contains also weight of fixed modification
                        self.weight -= AminoAcid::get(aa_one_letter_code).get_mono_mass();
                        self.weight += AminoAcid::get(*aa_replacement).get_mono_mass();
                        match self.aa_sequence.get_mut(idx) {
                            Some(item) => *item = *aa_replacement,
                            None => panic!("proteomic::models::decoys::new_decoy::NewDecoy.swap_amino_acids_to_hit_mass_tolerance(): expected char insted of None at self.aa_sequence.get_mut(idx)")
                        }
                        if let Some(ref modification) = fix_modifications_map.get(aa_replacement) {
                            match self.add_modification_at(idx, modification) {
                                Ok(_) => (),
                                Err(err) => panic!("proteomic::models::decoys::new_decoy::NewDecoy.swap_amino_acids_to_hit_mass_tolerance(): {}", err)
                            };
                        }
                        if self.hits_mass_tolerance() {
                            if self.create(conn) { created_counter += 1 }
                        }
                        created_counter += self.try_variable_modifications(conn, max_number_of_modifications, varibale_modification_map);
                        continue 'sequence; // continue with next amino acid in aa_sequence after swapping
                    }
                }
            }
            tries += 1;
            if tries == 1000 { break 'tries }
        }
        return created_counter;
    }

    fn try_variable_modifications(&mut self, conn: &postgres::Connection, max_number_of_modifications: u8, varibale_modification_map: &HashMap<char, Modification>) -> usize {
        let mut created_counter: usize = 0;
        let mut modification_positions: Vec<(usize, char)> = Vec::new();
        for (idx, one_letter_code) in self.aa_sequence.iter().enumerate() {
            if varibale_modification_map.contains_key(one_letter_code) {
                modification_positions.push((idx, *one_letter_code));
            }
        }
        for number_of_modifications in 1..=max_number_of_modifications {
            if number_of_modifications as usize <= modification_positions.len() {
                let n_choose_k = NChooseK::new(number_of_modifications as i32, modification_positions.clone());
                'combinations: for combination in n_choose_k.into_iter() {
                    self.remove_all_variable_modifications();
                    'positions: for position in combination {
                        if let Some(modification) = varibale_modification_map.get(&position.1) {
                            match self.set_variable_modification_at(position.0, modification) {
                                Ok(_) => continue 'positions,
                                Err(err) => match err {
                                    NewDecoyError::OutrangeMassTolerance => continue 'combinations,
                                    NewDecoyError::AlreadyFixModificationInPlace => continue 'positions,
                                    _ => panic!("proteomic::models::decoys::new_decoy::NewDecoy.try_variable_modifications(): {}", err)
                                }
                            }
                        }
                    }
                    if self.hits_mass_tolerance() {
                        if self.create(conn) { created_counter += 1 }
                    }
                }
            }
        }
        return created_counter;
    }

    fn create(&self, conn: &postgres::Connection) -> bool {
        let mut result = false;
        if self.hits_mass_tolerance() {
            let mut base_decoy: BaseDecoy = self.as_base_decoy();
            // check if decoy is peptide
            match Peptide::exists_where(&conn, "aa_sequence = $1", &[&self.get_aa_sequence()]) {
                Ok(_) => return false,
                Err(err) => match err {
                    QueryError::NoMatch => (),  // if no match is found continue with this decoy
                    _ => panic!("proteomic::utility::decoy_generator.generate_decoys() could not check if decoy is peptide: {}", err)
                }
            }
            match base_decoy.create(&conn) {
                Ok(query_ok) => match query_ok {
                    QueryOk::Created => {
                        if self.number_of_modifications == 0 { result = true; }     // base decoy should only be counted if no modifications are applied
                    },
                    QueryOk::AlreadyExists => (),
                    _ => panic!("proteomic::utility::decoy_generator::generate(): In fact not other QueryOk than QueryOk::Created and QueryOk::AlreadyExists are used in BaseDecoy.create(), so this panic shoud not be reached.")
                },
                Err(err) => println!("proteomic::utility::decoy_generator::generate() when second BaseDecoy.create(): {}", err)
            }
            if self.number_of_modifications > 0 {
                let mut modified_decoy = self.as_modified_decoy(&base_decoy);
                match modified_decoy.create(&conn) {
                    Ok(query_ok) => match query_ok {
                        QueryOk::Created => result = true,
                        QueryOk::AlreadyExists => (),
                        _ => panic!("proteomic::utility::decoy_generator::generate(): In fact not other QueryOk than QueryOk::Created and QueryOk::AlreadyExists are used in BaseDecoy.create(), so this panic shoud not be reached.")
                    },
                    Err(err) => println!("proteomic::utility::decoy_generator::generate() when second BaseDecoy.create(): {}", err)
                }
            }
        }
        return result;
    }
}

impl Decoy for NewDecoy {
    fn to_string(&self) -> String {
        return format!(
            "proteomic::modes::decoys::new_decoy::NewDecoy\n\taa_sequence => {}\n\tweight => {}\n\tc-terminus-modification => {}\n\tn-terminus-modification => {}\n\tmodifications => {}",
            self.get_aa_sequence(),
            mass::convert_mass_to_float(self.weight),
            match self.c_terminus_modification {
                Some(ref modification) => modification.get_accession(),
                None => "None"
            },
            match self.n_terminus_modification {
                Some(ref modification) => modification.get_accession(),
                None => "None"
            },
            self.modifications.iter().map(
                |modification_option| match modification_option {
                    Some(ref modification) => modification.get_accession().to_string(),
                    None => "None".to_string()
                }
            ).collect::<Vec<String>>().join(", ")
        );
    }

    fn get_header(&self) -> String {
        return ">decoy|NONE|None Randomized-decoy OS=None OX=None SV=None".to_owned();
    }

    fn get_aa_sequence(&self) -> String {
        return self.aa_sequence.iter().collect::<String>();
    }

    fn get_weight(&self) -> i64 {
        return self.weight;
    }

    fn get_length(&self) -> i32 {
        return self.aa_sequence.len() as i32;
    }

    fn get_c_terminus_amino_acid(&self) -> char {
        match self.aa_sequence.iter().last() {
            Some(amino_acids_one_letter_code) => *amino_acids_one_letter_code,
            None => '_'
        }
    }

    fn get_n_terminus_amino_acid(&self) -> char {
        match self.aa_sequence.iter().next() {
            Some(amino_acids_one_letter_code) => *amino_acids_one_letter_code,
            None => '_'
        }
    }

    fn get_amino_acid_at(&self, idx: usize) -> char {
        match self.aa_sequence.iter().nth(idx) {
            Some(amino_acids_one_letter_code) => *amino_acids_one_letter_code,
            None => '_'
        }
    }
}
