use proteomic::models::enzyms::transaction_summary::TransactionSummary;

pub struct DigestSummary {
    has_created_protein: bool,
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
            has_created_protein: false,
            number_of_created_peptides: 0,
            number_of_created_peptide_protein_associations: 0,
            number_of_processed_peptides: 0,
            number_of_processed_peptide_protein_associations: 0,
            unsolveable_errors_occured: false,
            log: Vec::new()
        }
    }

    pub fn has_created_protein(&self) -> bool {
        return self.has_created_protein;
    }

    pub fn set_has_created_protein(&mut self) {
        self.has_created_protein = true;
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

    pub fn increase_peptide_protein_association_counter(&mut self, association_created: bool) {
        if association_created { self.number_of_created_peptide_protein_associations += 1; }
        self.number_of_processed_peptide_protein_associations += 1;
    }

    pub fn increase_peptides_counter(&mut self, peptide_created: bool) {
        if peptide_created { self.number_of_created_peptides += 1; }
        self.number_of_processed_peptides += 1;
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

    pub fn merge_with_transaction_summary(&mut self, transaction_summary: &TransactionSummary) {
        self.number_of_created_peptides += transaction_summary.get_number_of_created_peptides();
        self.number_of_created_peptide_protein_associations += transaction_summary.get_number_of_created_peptide_protein_associations();
        self.number_of_processed_peptides += transaction_summary.get_number_of_processed_peptides();
        self.number_of_processed_peptide_protein_associations += transaction_summary.get_number_of_processed_peptide_protein_associations();
    }
}