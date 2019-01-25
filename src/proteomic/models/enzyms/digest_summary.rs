pub struct DigestSummary {
    created_protein: bool,
    number_of_created_peptides: usize,
    number_of_created_peptide_protein_associations: usize,
    number_of_processed_peptides: usize,
    number_of_processed_peptide_protein_associations: usize,
    unsolveable_errors_occured: bool,
    log: Vec<String>
}

impl DigestSummary {
    pub fn new() -> DigestSummary {
        return DigestSummary {
            created_protein: false,
            number_of_created_peptides: 0,
            number_of_created_peptide_protein_associations: 0,
            number_of_processed_peptides: 0,
            number_of_processed_peptide_protein_associations: 0,
            unsolveable_errors_occured: false,
            log: Vec::new()
        }
    }

    pub fn get_created_protein(&self) -> bool {
        return self.created_protein;
    }

    pub fn set_created_protein(&mut self) {
        self.created_protein = true;
    }

    pub fn get_number_of_created_peptides(&self) -> usize {
        return self.number_of_created_peptides;
    }

    pub fn get_number_of_created_peptide_protein_associations(&self) -> usize {
        return self.number_of_created_peptide_protein_associations;
    }

    pub fn get_number_of_processed_peptides(&self) -> usize {
        return self.number_of_processed_peptides;
    }

    pub fn get_number_of_processed_peptide_protein_associations(&self) -> usize {
        return self.number_of_processed_peptide_protein_associations;
    }

    pub fn increase_counter(&mut self, peptide_created: bool, association_created: bool) {
        if peptide_created { self.number_of_created_peptides += 1; }
        self.number_of_processed_peptides += 1;
        if association_created { self.number_of_created_peptide_protein_associations += 1; }
        self.number_of_processed_peptide_protein_associations += 1;
    }

    pub fn get_log(&self) -> String {
        return self.log.join("\n");
    }

    pub fn log_push(&mut self, message: &str) {
        self.log.push(message.to_owned());
    }

    pub fn get_number_of_log_messages(&self) ->  usize {
        return self.log.len();
    }

    pub fn get_unsolveable_errors_occured(&self) -> bool {
        return self.unsolveable_errors_occured;
    }

    pub fn set_unsolveable_errors_occured(&mut self) {
        self.unsolveable_errors_occured = true;
    }
}