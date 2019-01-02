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

    pub fn get_difference_between_modified_weight_and_upper_weight_limit(&self) -> i64 {
        return self.upper_weight_limit - self.modified_weight;
    }

    pub fn get_number_of_modifications(&self) -> i32 {
        return self.number_of_modifications;
    }

    pub fn hits_mass_tolerance(&self) -> bool {
        return ((self.lower_weight_limit <= self.weight) & (self.weight <= self.upper_weight_limit))
            | ((self.lower_weight_limit <= self.modified_weight) & (self.modified_weight <= self.upper_weight_limit));
    }

    pub fn difference_to_upper_weight_limit(&self) -> i64 {
        return self.upper_weight_limit - self.weight;
    }

    fn push_amino_acid(&mut self, amino_acid: &AminoAcid) {
        // after adding a new amino acid, the current has no free c-terminus anymore. so we have to clean the current c-terminus-modification
        if self.c_terminus_modification.is_some() {
            self.set_c_terminus_modification_to_none();
        }
        self.weight += amino_acid.get_mono_mass();
        self.modified_weight += amino_acid.get_mono_mass();
        self.aa_sequence.push(amino_acid.get_one_letter_code());
    }

    fn push_modification(&mut self, modification_option: &Option<&Modification>) {
        match modification_option {
            Some(modification) => {
                // if this is the first modification at all and the for the first amino acid
                if (self.modifications.len() == 0) & (self.get_length() == 1) {
                    match modification.get_position() {
                        // modification for n-terminus
                        ModificationPosition::NTerminus => {
                            self.modified_weight += modification.get_mono_mass();
                            self.n_terminus_modification = Some((*modification).clone());
                            self.number_of_modifications += 1;
                        }
                        // modification for anywhere
                        ModificationPosition::Anywhere => {
                            self.modified_weight += modification.get_mono_mass();
                            self.modifications.push(Some((*modification).clone()));
                            self.number_of_modifications += 1;
                        }
                        // else
                        _ => self.modifications.push(None)
                    }
                } else {
                    match modification.get_position() {
                        // modification for c-terminus
                        ModificationPosition::CTerminus => {
                            self.modified_weight += modification.get_mono_mass();
                            self.c_terminus_modification = Some((*modification).clone());
                            self.number_of_modifications += 1;
                        }
                        // modification for anywhere
                        ModificationPosition::Anywhere => {
                            self.modified_weight += modification.get_mono_mass();
                            self.modifications.push(Some((*modification).clone()));
                            self.number_of_modifications += 1;
                        }
                        // else
                        _ => self.modifications.push(None)
                    }
                }
            },
            None => self.modifications.push(None)
        }
    }

    pub fn push_amino_acid_with_modification(&mut self, amino_acid: &AminoAcid, modification_option: &Option<&Modification>) -> Result<(), NewDecoyError> {
        match modification_option {
            // check if modfication matches amino acid
            Some(modification) if amino_acid.get_one_letter_code() != modification.get_amino_acid_one_letter_code() => return Err(NewDecoyError::ModificationDoesNotMatchToAminoAcid),
            Some(_) => (),  // else
            None => (),
        }
        self.push_amino_acid(amino_acid);
        self.push_modification(modification_option);
        return Ok(());
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
            Some(modification) => {
                self.modified_weight -= modification.get_mono_mass();
                self.number_of_modifications -= 1;
            },
            None => ()
        };
        self.c_terminus_modification = None;
        return c_terminus_modification_option_clone;
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
