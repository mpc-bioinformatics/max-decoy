use proteomic::models::amino_acids::modification::{Modification, ModificationPosition};
use proteomic::models::amino_acids::amino_acid::AminoAcid;
use proteomic::models::mass;
use proteomic::models::mass::neutral_loss::NeutralLoss;
use proteomic::models::decoys::decoy::Decoy;
use proteomic::models::decoys::base_decoy::BaseDecoy;
use proteomic::models::decoys::modified_decoy::ModifiedDecoy;

pub enum NewDecoyError {
    ModificationDoesNotMatchToAminoAcid        // modification does not relate to amino acid
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

    pub fn get_modified_weight(&self) -> i64 {
        return self.modified_weight;
    }

    pub fn get_number_of_modifications(&self) -> i32 {
        return self.number_of_modifications;
    }

    pub fn hits_mass_tolerance(&self) -> bool {
        return (self.lower_weight_limit <= self.modified_weight) & (self.modified_weight <= self.upper_weight_limit);
    }

    pub fn get_distance_to_hit_mass_tolerance(&self) -> i64 {
        if self.modified_weight < self.lower_weight_limit {
            return self.lower_weight_limit - self.modified_weight;
        } else if self.modified_weight > self.lower_weight_limit {
            return self.modified_weight - self.upper_weight_limit;
        } else {
            return 0;
        }
    }

    pub fn push_amino_acid(&mut self, amino_acid: &AminoAcid) {
        self.weight += amino_acid.get_mono_mass();
        self.modified_weight += amino_acid.get_mono_mass();
        self.aa_sequence.push(amino_acid.get_one_letter_code());
        self.modifications.push(None);
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

    pub fn set_n_terminus_modification(&mut self, modification: &Modification) -> Result<(), NewDecoyError> {
        if modification.get_amino_acid_one_letter_code() != self.get_n_terminus_amino_acid() {
            return Err(NewDecoyError::ModificationDoesNotMatchToAminoAcid);
        }
        match &self.n_terminus_modification {
            Some(modification) => self.modified_weight -= modification.get_mono_mass(),
            None => ()
        }
        self.modified_weight += modification.get_mono_mass();
        self.n_terminus_modification = Some(modification.clone());
        self.number_of_modifications += 1;
        return Ok(());
    }

    pub fn set_c_terminus_modification(&mut self, modification: &Modification) -> Result<(), NewDecoyError> {
        if modification.get_amino_acid_one_letter_code() != self.get_c_terminus_amino_acid() {
            return Err(NewDecoyError::ModificationDoesNotMatchToAminoAcid);
        }
        match &self.c_terminus_modification {
            Some(modification) => self.modified_weight -= modification.get_mono_mass(),
            None => ()
        }
        self.modified_weight += modification.get_mono_mass();
        self.c_terminus_modification = Some(modification.clone());
        self.number_of_modifications += 1;
        return Ok(());
    }

    pub fn set_modification_at(&mut self, idx: usize, modification: &Modification) -> Result<(), NewDecoyError> {
        if modification.get_amino_acid_one_letter_code() != self.get_amino_acid_at(idx) {
            return Err(NewDecoyError::ModificationDoesNotMatchToAminoAcid);
        }
        match self.modifications.get_mut(idx) {
            Some(modification_option) => {
                match modification_option {
                    Some(modification) => self.modified_weight -= modification.get_mono_mass(),
                    None => ()
                }
                self.modified_weight += modification.get_mono_mass();
                *modification_option = Some(modification.clone());
                self.number_of_modifications += 1;
            }
            None => ()
        }
        return Ok(());
    }

    pub fn remove_all_modifications(&mut self) {
        match &self.n_terminus_modification {
            Some(modification) => self.modified_weight -= modification.get_mono_mass(),
            None => ()
        }
        self.n_terminus_modification = None;
        match &self.c_terminus_modification {
            Some(modification) => self.modified_weight -= modification.get_mono_mass(),
            None => ()
        }
        for modification_option in self.modifications.iter_mut() {
            match modification_option {
                Some(modification) => self.modified_weight -= modification.get_mono_mass(),
                None => ()
            }
            *modification_option = None;
        }
        self.number_of_modifications = 0;
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
