use std::collections::HashMap;

use proteomic::models::amino_acids::modification::{Modification, ModificationPosition};
use proteomic::models::amino_acids::amino_acid::AminoAcid;
use proteomic::models::mass;
use proteomic::models::mass::neutral_loss::NeutralLoss;
use proteomic::models::decoys::decoy::Decoy;
use proteomic::models::decoys::base_decoy::BaseDecoy;
use proteomic::models::decoys::modified_decoy::ModifiedDecoy;

pub enum NewDecoyError {
    ModificationDoesNotMatchToAminoAcid,       // modification does not relate to amino acid
    OutrangeMassTolerance,
    ModificationIsNotFix,
    ModificationIsNotVariable,
    AlreadyFixModificationInPlace
}

pub struct NewDecoy {
    aa_sequence: Vec<char>,
    modifications: Vec<Option<Modification>>,
    weight: i64,
    modified_weight: i64,
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
            modified_weight: water_loss.get_mono_mass(),
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

    pub fn get_modified_weight(&self) -> i64 {
        return self.modified_weight;
    }

    pub fn get_number_of_modifications(&self) -> i32 {
        return self.number_of_modifications;
    }

    pub fn hits_mass_tolerance(&self) -> bool {
        return (self.lower_weight_limit <= self.modified_weight) & (self.modified_weight <= self.upper_weight_limit);
    }

    /// returns positive value if modified_weight is less than lower_weight_limit
    /// returns negative value if modified_weight is higher than upper_weight_limit
    /// returns zero if modified_weight is between or equals lower_weight_limit and upper_weight_limit
    pub fn get_distance_to_mass_tolerance(&self) -> i64 {
        if self.modified_weight < self.lower_weight_limit {
            return self.lower_weight_limit - self.modified_weight;
        } else if self.modified_weight > self.upper_weight_limit {
            return self.upper_weight_limit - self.modified_weight;
        } else {
            return 0;
        }
    }

    fn push_amino_acid(&mut self, amino_acid: &AminoAcid) {
        self.weight += amino_acid.get_mono_mass();
        self.modified_weight += amino_acid.get_mono_mass();
        self.aa_sequence.push(amino_acid.get_one_letter_code());
    }

    fn undo_push_amino_acid(&mut self, amino_acid: &AminoAcid) {
        self.weight -= amino_acid.get_mono_mass();
        self.modified_weight -= amino_acid.get_mono_mass();
        let last_idx = self.aa_sequence.len() - 1;
        self.aa_sequence.remove(last_idx);
    }

    /// returns the old c_terminus_modification
    fn push_modification(&mut self, modification_option: &Option<Modification>) -> Option<Modification> {
        let old_c_terminus_modification = self.remove_c_terminus_modification();        // after pushing a new amino acid, the old-last one has no c-terminus any more. so remove the c-terminus modification
        if let Some(modification) = modification_option {
            match modification.get_position() {
                ModificationPosition::Anywhere => {
                    self.modified_weight += modification.get_mono_mass();
                    self.number_of_modifications += 1;
                    self.modifications.push(Some(modification.clone()));
                },
                ModificationPosition::CTerminus => {
                    self.modified_weight += modification.get_mono_mass();
                    self.number_of_modifications += 1;
                    self.c_terminus_modification = Some(modification.clone());
                },
                // only if this is the first amino acid and modification
                ModificationPosition::NTerminus if self.aa_sequence.len() == 1 => {
                    self.modified_weight += modification.get_mono_mass();
                    self.number_of_modifications += 1;
                    self.n_terminus_modification = Some(modification.clone());
                },
                _ => ()
            }
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
                    self.modified_weight -= modification.get_mono_mass();
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
                    self.modified_weight += modification.get_mono_mass();
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
            self.weight
        )
    }

    pub fn as_modified_decoy(&self, base_decoy: &BaseDecoy) -> ModifiedDecoy {
        return ModifiedDecoy::new(
            base_decoy,
            self.c_terminus_modification.clone(),
            self.n_terminus_modification.clone(),
            &self.modifications,
            self.modified_weight
        );
    }

    fn remove_c_terminus_modification(&mut self) -> Option<Modification> {
        match self.c_terminus_modification {
            Some(ref modification) => {
                self.modified_weight -= modification.get_mono_mass();
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
                self.modified_weight -= modification.get_mono_mass();
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

    pub fn set_variable_c_terminus_modification(&mut self, modification: &Modification) -> Result<(), NewDecoyError> {
        if !modification.is_fix() {
            return Err(NewDecoyError::ModificationIsNotVariable);
        }
        if modification.get_amino_acid_one_letter_code() != self.get_c_terminus_amino_acid() {
            return Err(NewDecoyError::ModificationDoesNotMatchToAminoAcid);
        }
        if self.c_terminus_modification.is_some() {
            return Err(NewDecoyError::AlreadyFixModificationInPlace);
        }
        self.modified_weight += modification.get_mono_mass();
        self.number_of_modifications += 1;
        self.c_terminus_modification = Some(modification.clone());
        if self.get_distance_to_mass_tolerance() < 0 {
            return Err(NewDecoyError::OutrangeMassTolerance);
        }
        return Ok(());
    }

    pub fn set_variable_n_terminus_modification(&mut self, modification: &Modification) -> Result<(), NewDecoyError> {
        if !modification.is_fix() {
            return Err(NewDecoyError::ModificationIsNotVariable);
        }
        if modification.get_amino_acid_one_letter_code() != self.get_n_terminus_amino_acid() {
            return Err(NewDecoyError::ModificationDoesNotMatchToAminoAcid);
        }
        if self.n_terminus_modification.is_some() {
            return Err(NewDecoyError::AlreadyFixModificationInPlace);
        }
        self.modified_weight += modification.get_mono_mass();
        self.number_of_modifications += 1;
        self.n_terminus_modification = Some(modification.clone());
        if self.get_distance_to_mass_tolerance() < 0 {
            return Err(NewDecoyError::OutrangeMassTolerance);
        }
        return Ok(());
    }

    pub fn set_variable_modification_at(&mut self, idx: usize, modification: &Modification) -> Result<(), NewDecoyError> {
        if !modification.is_fix() {
            return Err(NewDecoyError::ModificationIsNotVariable);
        }
        if modification.get_amino_acid_one_letter_code() != self.get_amino_acid_at(idx) {
            return Err(NewDecoyError::ModificationDoesNotMatchToAminoAcid);
        }
        match self.modifications.get(idx) {
            Some(modification_option) => match modification_option {
                Some(modification) => return Err(NewDecoyError::AlreadyFixModificationInPlace),
                None => ()
            },
            None => ()
        }
        self.modified_weight += modification.get_mono_mass();
        self.number_of_modifications += 1;
        self.modifications[idx] = Some(modification.clone());
        if self.get_distance_to_mass_tolerance() < 0 {
            return Err(NewDecoyError::OutrangeMassTolerance);
        }
        return Ok(());
    }

    pub fn remove_all_variable_modifications(&mut self) {
        match self.n_terminus_modification {
            Some(ref mut modification) if !modification.is_fix() => {
                self.modified_weight -= modification.get_mono_mass();
                self.number_of_modifications -= 1;
            },
            _ => ()
        }
        self.n_terminus_modification = None;

        match self.c_terminus_modification {
            Some(ref mut modification) if !modification.is_fix() => {
                self.modified_weight -= modification.get_mono_mass();
                self.number_of_modifications -= 1;
            },
            _ => ()
        }
        self.c_terminus_modification = None;

        for modification_option in self.modifications.iter_mut() {
            match modification_option {
                Some(modification) if !modification.is_fix() => {
                    self.modified_weight -= modification.get_mono_mass();
                    self.number_of_modifications -= 1;
                },
                _ => ()
            }
            *modification_option = None;
        }
    }
}

impl Decoy for NewDecoy {
    fn to_string(&self) -> String {
        return format!(
            "proteomic::modes::decoys::new_decoy::NewDecoy\n\taa_sequence => {}\n\tweight => {}",
            self.get_aa_sequence(),
            mass::convert_mass_to_float(self.weight)
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
