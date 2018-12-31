use proteomic::models::amino_acids::modification::{Modification, ModificationPosition};
use proteomic::models::amino_acids::amino_acid::AminoAcid;
use proteomic::models::mass;
use proteomic::models::mass::neutral_loss::NeutralLoss;
use proteomic::models::persistable::Persistable;
use proteomic::models::decoys::decoy::Decoy;
use proteomic::models::decoys::base_decoy::BaseDecoy;
use proteomic::models::decoys::modified_decoy::ModifiedDecoy;

pub struct NewDecoy {
    aa_sequence: Vec<char>,
    modifications: Vec<Option<Modification>>,
    weight: i64,
    modified_weight: i64,
    lower_weight_limit: i64,
    upper_weight_limit: i64,
    n_terminus_modification: Option<Modification>,
    c_terminus_modification: Option<Modification>
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
            c_terminus_modification: None
        }
    }

    pub fn get_modified_weight(&self) -> i64 {
        return self.modified_weight;
    }

    pub fn hits_mass_tolerance(&self) -> bool {
        return ((self.lower_weight_limit <= self.weight) & (self.weight <= self.upper_weight_limit))
            | ((self.lower_weight_limit <= self.modified_weight) & (self.modified_weight <= self.upper_weight_limit));
    }

    pub fn difference_to_upper_weight_limit(&self) -> i64 {
        return self.upper_weight_limit - self.weight;
    }

    fn add_amino_acid(&mut self, amino_acid: &AminoAcid) -> bool {
        let temp_weight: i64 = self.weight + amino_acid.get_mono_mass();
        if self.is_amino_acid_fitting(amino_acid) {
            self.weight += amino_acid.get_mono_mass();
            self.aa_sequence.push(amino_acid.get_one_letter_code());
            self.modifications.push(None);
            return true;
        } else {
            return false;
        }
    }

    fn is_amino_acid_fitting(&self, amino_acid: &AminoAcid) -> bool {
        let temp_weight: i64 = self.weight + amino_acid.get_mono_mass();
        return (self.lower_weight_limit <= temp_weight) & (temp_weight <= self.upper_weight_limit)
    }

    fn is_modification_fitting(&self, modification: &Modification) -> bool {
        let temp_weight: i64 = self.modified_weight + modification.get_mono_mass();
        return (self.lower_weight_limit <= temp_weight) & (temp_weight <= self.upper_weight_limit);
    }

    // removes the old modification and returns it
    fn set_c_terminus_modification_to_none(&mut self) -> Option<Modification> {
        let c_terminus_modification_option_clone: Option<Modification> = self.c_terminus_modification.clone();
        match &self.c_terminus_modification {
            Some(modification) => self.modified_weight -= modification.get_mono_mass(),
            None => ()
        };
        self.c_terminus_modification = None;
        return c_terminus_modification_option_clone;
    }

    fn set_n_terminus_modification_to_none(&mut self)  -> Option<Modification> {
        let n_terminus_modification_option_clone: Option<Modification> = self.c_terminus_modification.clone();
        match &self.n_terminus_modification {
            Some(modification) => self.modified_weight -= modification.get_mono_mass(),
            None => ()
        };
        self.n_terminus_modification = None;
        return n_terminus_modification_option_clone;
    }

    fn set_modification_to_none_at(&mut self, idx: usize)  -> Option<Modification> {
        let modification_option_clone: Option<Modification> = self.modifications.get(idx).unwrap().clone();
        match &mut self.modifications.get_mut(idx) {
            Some(ref mut modification_option) => {
                match modification_option {
                    Some(modification) => {
                        self.modified_weight -= modification.get_mono_mass()

                    },
                    None => ()
                }
                **modification_option = None;
            }
            None => ()
        };
        return modification_option_clone;
    }

    pub fn set_c_terminus_modification(&mut self, new_c_terminus_modification_option: Option<Modification>) -> Result<(), String> {
        match new_c_terminus_modification_option {
            Some(modification) => {
                if modification.get_position() != ModificationPosition::CTerminus {
                    return Err("Modification is not for C-Terminus".to_owned());
                }
                if modification.get_amino_acid_one_letter_code() == self.get_c_terminus_amino_acid() {
                    return Err(format!("Modification is not for the amino acid on C-Terminus: {}", self.get_c_terminus_amino_acid()))
                }
                let old_modification: Option<Modification> = self.set_c_terminus_modification_to_none();
                if self.is_modification_fitting(&modification) {
                    self.modified_weight += modification.get_mono_mass();
                    self.c_terminus_modification = Some(modification);
                } else {
                    match &old_modification {
                        Some(old_mod) => self.modified_weight += old_mod.get_mono_mass(),
                        None => ()
                    }
                    self.c_terminus_modification = old_modification;
                }
                return Ok(());
            }
            None => {
                self.set_c_terminus_modification_to_none();
                return Ok(());
            }
        }
    }

    pub fn set_n_terminus_modification(&mut self, new_n_terminus_modification: Option<Modification>) -> Result<(), String> {
        match new_n_terminus_modification {
            Some(modification) => {
                if modification.get_position() != ModificationPosition::NTerminus {
                    return Err("Modification is not for N-Terminus".to_owned());
                }
                if modification.get_amino_acid_one_letter_code() == self.get_n_terminus_amino_acid() {
                    return Err(format!("Modification is not for the amino acid on N-Terminus: {}", self.get_n_terminus_amino_acid()))
                }
                let old_modification: Option<Modification> = self.set_n_terminus_modification_to_none();
                if self.is_modification_fitting(&modification) {
                    self.modified_weight += modification.get_mono_mass();
                    self.n_terminus_modification = Some(modification);
                } else {
                    match &old_modification {
                        Some(old_mod) => self.modified_weight += old_mod.get_mono_mass(),
                        None => ()
                    }
                    self.n_terminus_modification = old_modification;
                }
                return Ok(());
            }
            None => {
                self.set_n_terminus_modification_to_none();
                return Ok(());
            }
        }
    }

    pub fn set_modification_at(&mut self, idx: usize, modification_option: Option<Modification>) -> Result<(), String> {
        match modification_option {
            Some(modification) => {
                if modification.get_position() != ModificationPosition::Anywhere {
                    return Err("Modification is not for Anywhere".to_owned());
                }
                if modification.get_amino_acid_one_letter_code() != self.get_amino_acid_at(idx) {
                    return Err(format!("Modification is not for the amino acid on index {}: {}", idx, self.get_n_terminus_amino_acid()))
                }
                let old_modification: Option<Modification> = self.set_modification_to_none_at(idx);
                if self.is_modification_fitting(&modification) {
                    self.modified_weight += modification.get_mono_mass();
                    self.modifications[idx] = Some(modification);
                } else {
                    match &old_modification {
                        Some(old_mod) => self.modified_weight += old_mod.get_mono_mass(),
                        None => ()
                    }
                    self.modifications[idx] = old_modification;
                }
                return Ok(());
            }
            None => {
                self.modifications[idx] = None;
                return Ok(());
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

    pub fn as_modified_decoy(&self, base_decoy: &BaseDecoy) -> Result<ModifiedDecoy, String> {
        if base_decoy.get_aa_sequence() != self.get_aa_sequence() {
            return Err("given BaseDecoy has not the same amino acid sequence as the NewDecoy".to_owned());
        }
        if base_decoy.is_persisted() {
            return Err("given BaseDecoy is not persisted".to_owned());
        }
        return Ok(
            ModifiedDecoy::new(
                base_decoy,
                self.c_terminus_modification.clone(),
                self.n_terminus_modification.clone(), 
                &self.modifications,
                self.modified_weight
            )
        );
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
        return ">decoy|NONE|None Randomized-decoy OS=None OX=None SV=None ModResUniMod={}".to_owned();
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
