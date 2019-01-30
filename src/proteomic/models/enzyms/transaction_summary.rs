pub struct TransactionSummary {
    number_of_created_peptides: usize,
    number_of_created_peptide_protein_associations: usize,
    number_of_processed_peptides: usize,
    number_of_processed_peptide_protein_associations: usize,
}

impl TransactionSummary {
    pub fn new() -> TransactionSummary {
        return TransactionSummary {
            number_of_created_peptides: 0,
            number_of_created_peptide_protein_associations: 0,
            number_of_processed_peptides: 0,
            number_of_processed_peptide_protein_associations: 0
        }
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
}