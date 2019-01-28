use std::collections::HashMap;
use std::fmt;

use rand::{thread_rng, Rng};
use rand::seq::SliceRandom;

use proteomic::models::amino_acids::modification::{Modification, ModificationPosition};
use proteomic::models::amino_acids::amino_acid::AminoAcid;
use proteomic::models::mass;
use proteomic::models::mass::neutral_loss::NeutralLoss;
use proteomic::models::peptides::peptide_interface::PeptideInterface;
use proteomic::models::decoys::base_decoy::BaseDecoy;
use proteomic::models::decoys::modified_decoy::ModifiedDecoy;
use proteomic::utility::combinations::n_choose_k::NChooseK;
use proteomic::models::persistable::{Persistable, QueryOk, QueryError};
use proteomic::models::peptides::peptide::Peptide;
use proteomic::models::amino_acids::amino_acid::AMINO_ACIDS_FOR_DECOY_GENERATION;

const DECOY_HEADER: &'static str = ">DECOY_";

pub enum PushAminoAcidOk {
    LessThenMassTolerance,
    GreterThenMassTolerance,
    HitsMassTolerance
}

pub enum PushAminoAcidError {
    ModificationDoesNotMatchToAminoAcid,
    ModificationIsNotFix
}

impl PushAminoAcidError {
    pub fn to_string(&self) -> String {
        match self {
            PushAminoAcidError::ModificationDoesNotMatchToAminoAcid => format!("PushAminoAcidError::ModificationDoesNotMatchToAminoAcid"),
            PushAminoAcidError::ModificationIsNotFix => format!("PushAminoAcidError::ModificationIsNotFix")
        }
    }
}

impl fmt::Display for PushAminoAcidError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        return write!(f, "{}", self.to_string());
    }
}


pub enum ModifiedPeptideError {
    ModificationDoesNotMatchToAminoAcid,       // modification does not relate to amino acid
    OutrangeMassTolerance,
    ModificationIsNotFix,
    ModificationIsNotVariable,
    AlreadyFixModificationInPlace
}

impl ModifiedPeptideError {
    pub fn to_string(&self) -> String {
        match self {
            ModifiedPeptideError::ModificationDoesNotMatchToAminoAcid => format!("ModifiedPeptideError::ModificationDoesNotMatchToAminoAcid"),
            ModifiedPeptideError::OutrangeMassTolerance => format!("ModifiedPeptideError::OutrangeMassTolerance"),
            ModifiedPeptideError::ModificationIsNotFix => format!("ModifiedPeptideError::ModificationIsNotFix"),
            ModifiedPeptideError::ModificationIsNotVariable => format!("ModifiedPeptideError::ModificationIsNotVariable"),
            ModifiedPeptideError::AlreadyFixModificationInPlace => format!("ModifiedPeptideError::AlreadyFixModificationInPlace")
        }
    }
}

impl fmt::Display for ModifiedPeptideError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        return write!(f, "{}", self.to_string());
    }
}



pub struct ModifiedPeptide {
    header: String,
    aa_sequence: Vec<char>,
    modifications: Vec<Option<Modification>>,
    weight: i64,
    precursor_mass: i64,
    lower_weight_limit: i64,
    upper_weight_limit: i64,
    n_terminus_modification: Option<Modification>,
    c_terminus_modification: Option<Modification>,
    number_of_modifications: i32
}

impl ModifiedPeptide {
    pub fn new(header: &str, precursor_mass: i64, lower_weight_limit: i64, upper_weight_limit: i64) -> Self {
        let water_loss: NeutralLoss = NeutralLoss::get("H2O");
        return Self {
            header: header.to_owned(),
            aa_sequence: Vec::new(),
            modifications: Vec::new(),
            weight: water_loss.get_mono_mass(),
            precursor_mass: precursor_mass,
            lower_weight_limit: lower_weight_limit,
            upper_weight_limit: upper_weight_limit,
            n_terminus_modification: None,
            c_terminus_modification: None,
            number_of_modifications: 0
        }
    }

    /// Creates a new ModifiedPepitide, with header DECOY_HEADER
    pub fn new_decoy(precursor_mass: i64, lower_weight_limit: i64, upper_weight_limit: i64) -> Self {
        return Self::new(DECOY_HEADER, precursor_mass, lower_weight_limit, upper_weight_limit);
    }

    /// Creates a new ModifiedPepitide from an existing sequence,
    /// with header DECOY_HEADER
    /// and respect to fix modifications
    pub fn decoy_from_string(aa_sequence: &str, precursor_mass: i64, lower_weight_limit: i64, upper_weight_limit: i64, fix_modifications: &HashMap<char, Modification>) -> Self {
        return Self::from_string(DECOY_HEADER, aa_sequence, precursor_mass, lower_weight_limit, upper_weight_limit, fix_modifications);
    }

    /// Creates a new ModifiedPepitide from an existing sequence and header
    /// with respect to fix modifications
    pub fn from_string(header: &str, aa_sequence: &str, precursor_mass: i64, lower_weight_limit: i64, upper_weight_limit: i64, fix_modifications: &HashMap<char, Modification>) -> Self {
        let mut new_modified_peptide = Self::new(header, precursor_mass, lower_weight_limit, upper_weight_limit);
        for amino_acids_one_letter_code in aa_sequence.chars() {
            let amino_acid = AminoAcid::get(amino_acids_one_letter_code);
            let modification_option = match fix_modifications.get(&amino_acids_one_letter_code) {
                Some(ref modification) => Some((*modification).clone()),
                None => None
            };
            match new_modified_peptide.push_amino_acid_and_fix_modification(&amino_acid, &modification_option) {
                Ok(_) => (),
                Err(err) => match err {
                    PushAminoAcidError::ModificationIsNotFix => panic!("proteomic::models::decoys::decoy::modified_peptide::ModifiedPeptide::from_string(): variable modification is passed to push_amino_acid_and_fix_modification(), which actually should not happen here."),
                    _ => continue
                }
            }
        }
        return new_modified_peptide;
    }

    /// Creates a new ModifiedPepitide from a Peptide
    /// with respect to fix modifications
    pub fn from_peptide(peptide: &Peptide, precursor_mass: i64, lower_weight_limit: i64, upper_weight_limit: i64, fix_modifications: &HashMap<char, Modification>) -> Self {
        return Self::from_string(">PEPTIDE ", peptide.get_aa_sequence(), precursor_mass, lower_weight_limit, upper_weight_limit, fix_modifications);
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

    /// Returns the distance to weight_limit
    pub fn get_distance_to_precursor_mass(&self) -> i64 {
        return (self.precursor_mass - self.weight).abs();
    }

    /// Pushs a new amino acid to sequence
    fn push_amino_acid(&mut self, amino_acid: &AminoAcid) {
        self.weight += amino_acid.get_mono_mass();
        self.aa_sequence.push(amino_acid.get_one_letter_code());
    }

    /// Removes the last amino acid.
    fn undo_push_amino_acid(&mut self, amino_acid: &AminoAcid) {
        self.weight -= amino_acid.get_mono_mass();
        let last_idx = self.aa_sequence.len() - 1;
        self.aa_sequence.remove(last_idx);
    }

    /// Pushs a new modification to the on od `self.modifications`.
    /// Returns the old c_terminus_modification in case you have to undo the push.
    /// ATTENTION: Use `self.push_amino_acid()` before, to make sure the modification fits the last amino acid
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

    /// Removes the last modification from `self.modification` and applys the old_c_terminus_modificiation to the 'new' end
    /// ATTENTION: Use `self.undo_push_amino_acid()` before, to make sure the last amino acid is removed
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
            self.get_aa_sequence().as_str()
        )
    }

    pub fn as_modified_decoy(&self, base_decoy: &BaseDecoy) -> ModifiedDecoy {
        let mut modifications: Vec<Modification> = Vec::new();
        if let Some(ref modification) = self.n_terminus_modification {
            modifications.push(modification.clone());
        }
        for modification_option in self.modifications.iter() {
            if let Some(ref modification) = modification_option {
                modifications.push(modification.clone());
            }
        }
        if let Some(ref modification) = self.c_terminus_modification {
            modifications.push(modification.clone());
        }
        return ModifiedDecoy::new(
            base_decoy,
            &modifications,
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

    /// Adds amino acid with modification.
    ///
    /// # Arguments
    ///
    /// * `amino_acid` - Amino acid
    /// * `modification_option` - Some(Modification )
    pub fn push_amino_acid_and_fix_modification(&mut self, amino_acid: &AminoAcid, modification_option: &Option<Modification>) -> Result<PushAminoAcidOk, PushAminoAcidError> {
        match modification_option {
            Some(modification) => {
                if !modification.is_fix() {
                    return Err(PushAminoAcidError::ModificationIsNotFix);
                }
                if modification.get_amino_acid_one_letter_code() != amino_acid.get_one_letter_code() {
                    return Err(PushAminoAcidError::ModificationDoesNotMatchToAminoAcid);
                }
            },
            None => ()
        }
        self.push_amino_acid(amino_acid);
        self.push_modification(modification_option);
        if self.weight < self.upper_weight_limit {
            return Ok(PushAminoAcidOk::LessThenMassTolerance);
        } else if self.weight > self.upper_weight_limit {
            return Ok(PushAminoAcidOk::GreterThenMassTolerance);
        } else {
            return Ok(PushAminoAcidOk::HitsMassTolerance);
        }
    }

    fn set_variable_c_terminus_modification(&mut self, modification: &Modification) -> Result<(), ModifiedPeptideError> {
        if modification.is_fix() {
            return Err(ModifiedPeptideError::ModificationIsNotVariable);
        }
        if modification.get_amino_acid_one_letter_code() != self.get_c_terminus_amino_acid() {
            return Err(ModifiedPeptideError::ModificationDoesNotMatchToAminoAcid);
        }
        if self.c_terminus_modification.is_some() {
            return Err(ModifiedPeptideError::AlreadyFixModificationInPlace);
        }
        self.weight += modification.get_mono_mass();
        self.number_of_modifications += 1;
        self.c_terminus_modification = Some(modification.clone());
        // if self.get_distance_to_mass_tolerance() < 0 {
        //     return Err(ModifiedPeptideError::OutrangeMassTolerance);
        // }
        return Ok(());
    }

    fn set_variable_n_terminus_modification(&mut self, modification: &Modification) -> Result<(), ModifiedPeptideError> {
        if modification.is_fix() {
            return Err(ModifiedPeptideError::ModificationIsNotVariable);
        }
        if modification.get_amino_acid_one_letter_code() != self.get_n_terminus_amino_acid() {
            return Err(ModifiedPeptideError::ModificationDoesNotMatchToAminoAcid);
        }
        if self.n_terminus_modification.is_some() {
            return Err(ModifiedPeptideError::AlreadyFixModificationInPlace);
        }
        self.weight += modification.get_mono_mass();
        self.number_of_modifications += 1;
        self.n_terminus_modification = Some(modification.clone());
        // if self.get_distance_to_mass_tolerance() < 0 {
        //     return Err(ModifiedPeptideError::OutrangeMassTolerance);
        // }
        return Ok(());
    }

    pub fn set_variable_modification_at(&mut self, idx: usize, modification: &Modification) -> Result<(), ModifiedPeptideError> {
        if modification.is_fix() {
            return Err(ModifiedPeptideError::ModificationIsNotVariable);
        }
        if modification.get_amino_acid_one_letter_code() != self.get_amino_acid_at(idx) {
            return Err(ModifiedPeptideError::ModificationDoesNotMatchToAminoAcid);
        }
        if (idx == 0) & (modification.get_position() == ModificationPosition::NTerminus) {
            return self.set_variable_n_terminus_modification(modification);
        } else if (idx == self.aa_sequence.len() - 1) & (modification.get_position() == ModificationPosition::CTerminus) {
            return self.set_variable_c_terminus_modification(modification);
        } else if modification.get_position() == ModificationPosition::Anywhere {
            match self.modifications.get_mut(idx) {
                Some(modification_option) => match modification_option {
                    Some(_) => return Err(ModifiedPeptideError::AlreadyFixModificationInPlace),
                    None => *modification_option = Some(modification.clone())
                },
                None => ()
            }
            self.weight += modification.get_mono_mass();
            self.number_of_modifications += 1;
            // if self.get_distance_to_mass_tolerance() < 0 {
            //     return Err(ModifiedPeptideError::OutrangeMassTolerance);
            // }
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

    fn add_modification_at(&mut self, idx: usize, modification: &Modification) -> Result<(), ModifiedPeptideError> {
        match self.aa_sequence.get(idx) {
            Some(one_letter_code) => {
                if *one_letter_code != modification.get_amino_acid_one_letter_code() {
                    return Err(ModifiedPeptideError::ModificationDoesNotMatchToAminoAcid);
                }
            },
            None => return Err(ModifiedPeptideError::ModificationDoesNotMatchToAminoAcid)
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
                None => panic!("proteomic::models::decoys::new_decoy::ModifiedPeptide.add_modification_at(): expect Some(Modification) instead of None")
            }
        }
        return Ok(());
    }

    pub fn swap_amino_acids_to_hit_mass_tolerance(&mut self, conn: &postgres::Connection, amino_acid_substitute_map: &HashMap<char, HashMap<char, i64>>, fix_modifications_map: &HashMap<char, Modification>, max_number_of_modifications: u8, varibale_modification_map: &HashMap<char, Modification>) -> bool {
        // try 100 times
        'tries: for _ in 0..100 {
            'sequence: for idx in 0..self.aa_sequence.len() {
                let aa_one_letter_code = self.get_amino_acid_at(idx);
                if let Some(ref swaps) = amino_acid_substitute_map.get(&aa_one_letter_code) {
                    let mut do_swap = false;
                    let mut best_swap: (i64, char) = (self.get_distance_to_precursor_mass(), aa_one_letter_code);
                    'swaps: for (amino_acid_swap, weight_change) in swaps.iter() {
                        self.weight += weight_change;   // temporarily increase weight
                        let distance = self.get_distance_to_precursor_mass();
                        if distance < best_swap.0 {
                            best_swap = (distance, *amino_acid_swap);
                            do_swap = true;

                        }
                        self.weight -= weight_change; // do not forget to reduce weight
                    }
                    if do_swap {
                        self.remove_modification_at(idx);
                        self.weight -= AminoAcid::get(aa_one_letter_code).get_mono_mass();
                        self.weight += AminoAcid::get(best_swap.1).get_mono_mass();
                        match self.aa_sequence.get_mut(idx) {
                            Some(item) => *item = best_swap.1,
                            None => panic!("proteomic::models::decoys::new_decoy::ModifiedPeptide.swap_amino_acids_to_hit_mass_tolerance(): expected char insted of None at self.aa_sequence.get_mut(idx)")
                        }
                        if let Some(ref modification) = fix_modifications_map.get(&best_swap.1) {
                            match self.add_modification_at(idx, modification) {
                                Ok(_) => (),
                                Err(err) => panic!("proteomic::models::decoys::new_decoy::ModifiedPeptide.swap_amino_acids_to_hit_mass_tolerance(): {}", err)
                            };
                        }
                        if self.hits_mass_tolerance() { return true; }
                        if self.try_variable_modifications(conn, max_number_of_modifications, varibale_modification_map) { return true; }
                    }
                }
            }
            // swap one random amino acid because gthe current sequence is at it' minimum
            let mut rng = thread_rng();
            let idx_to_swap: usize = rng.gen_range(0, self.aa_sequence.len());
            let aa_one_letter_code = self.get_amino_acid_at(idx_to_swap);
            let random_replacement = *AMINO_ACIDS_FOR_DECOY_GENERATION.choose(&mut rng).unwrap();
            self.remove_modification_at(idx_to_swap);
            self.weight -= AminoAcid::get(aa_one_letter_code).get_mono_mass();
            self.weight += AminoAcid::get(random_replacement).get_mono_mass();
            match self.aa_sequence.get_mut(idx_to_swap) {
                Some(item) => *item = random_replacement,
                None => panic!("proteomic::models::decoys::new_decoy::ModifiedPeptide.swap_amino_acids_to_hit_mass_tolerance(): expected char insted of None at self.aa_sequence.get_mut(idx_to_swap)")
            }
            if let Some(ref modification) = fix_modifications_map.get(&random_replacement) {
                match self.add_modification_at(idx_to_swap, modification) {
                    Ok(_) => (),
                    Err(err) => panic!("proteomic::models::decoys::new_decoy::ModifiedPeptide.swap_amino_acids_to_hit_mass_tolerance(): {}", err)
                };
            }
        }
        return false;
    }

    pub fn try_variable_modifications(&mut self, conn: &postgres::Connection, max_number_of_modifications: u8, varibale_modification_map: &HashMap<char, Modification>) -> bool {
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
                                    ModifiedPeptideError::OutrangeMassTolerance => continue 'combinations,
                                    ModifiedPeptideError::AlreadyFixModificationInPlace => continue 'positions,
                                    _ => panic!("proteomic::models::decoys::new_decoy::ModifiedPeptide.try_variable_modifications(): {}", err)
                                }
                            }
                        }
                    }
                    if self.hits_mass_tolerance() { return true; }
                }
            }
        }
        return false;
    }

    pub fn create(&self, conn: &postgres::Connection) -> bool {
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

impl PeptideInterface for ModifiedPeptide {
    fn to_string(&self) -> String {
        return format!(
            "proteomic::modes::decoys::new_decoy::ModifiedPeptide\n\taa_sequence => {}\n\tweight => {}\n\tc-terminus-modification => {}\n\tn-terminus-modification => {}\n\tmodifications => {}",
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
        return self.header.clone();
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
