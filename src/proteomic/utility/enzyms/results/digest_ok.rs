pub struct DigestOk {
    processed_peptides: usize,
    commited_peptides: usize,
    processed_peptide_protein_associations: usize,
    commited_peptide_protein_associations: usize,
    log_message: String
}

impl DigestOk {
    pub fn new(processed_peptides: usize, commited_peptides: usize, processed_peptide_protein_associations: usize, commited_peptide_protein_associations: usize, log_message: &str) -> DigestOk {
        return DigestOk {
            processed_peptides: processed_peptides,
            commited_peptides: commited_peptides,
            processed_peptide_protein_associations: processed_peptide_protein_associations,
            commited_peptide_protein_associations: commited_peptide_protein_associations,
            log_message: log_message.to_owned()
        }
    }

    pub fn get_processed_peptides(&self) -> usize {
        return self.processed_peptides;
    }

    pub fn get_commited_peptides(&self) -> usize {
        return self.commited_peptides;
    }

    pub fn get_processed_peptide_protein_associations(&self) -> usize {
        return self.processed_peptide_protein_associations;
    }

    pub fn get_commited_peptide_protein_associations(&self) -> usize {
        return self.commited_peptide_protein_associations;
    }

    pub fn get_log_message(&self) -> &String {
        return &self.log_message;
    }

}