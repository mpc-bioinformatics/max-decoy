/*
 * represents the position of a peptide in protein
 */

pub struct PeptidePosition {
    protein_id: i64,
    peptide_id: i64,
    position: usize
}

impl PeptidePosition {
    pub fn new(protein_id: i64, peptide_id: i64,  position: usize) -> PeptidePosition {
        return PeptidePosition {
            protein_id: protein_id,
            peptide_id: peptide_id,
            position: position
        }
    }
}