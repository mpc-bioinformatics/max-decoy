extern crate postgres;

use self::postgres::Connection;
use self::postgres::rows::Rows;
use self::postgres::Result;
use self::postgres::stmt::Statement;

use proteomic::models::peptide::Peptide;
use proteomic::models::protein::Protein;
use proteomic::models::persistable::Persistable;

pub struct PeptideProteinAssociation {
    peptide_id: i32,
    protein_id: i32
}

impl PeptideProteinAssociation {
    pub fn new(peptide: &Peptide, protein: &Protein) -> PeptideProteinAssociation {
        return PeptideProteinAssociation {
            peptide_id: peptide.get_id(),
            protein_id: protein.get_id(),
        }
    }

    pub fn get_peptide_id(&self) -> i32 {
        return self.peptide_id;
    }

    pub fn get_protein_id(&self) -> i32 {
        return self.protein_id;
    }

    pub fn create(&mut self, conn: &Connection) {
        self.execute_insert_query(&conn);
    }

    fn get_insert_statement() -> &'static str {
        return "INSERT INTO peptides_proteins (peptide_id, protein_id) VALUES ($1, $2) ON CONFLICT DO NOTHING";
    }

    fn execute_insert_query(&self, connection: &Connection) -> Result<Rows> {
        return connection.query(
            PeptideProteinAssociation::get_insert_statement(),
            &[&self.peptide_id, &self.protein_id]
        );
    }

    fn execute_insert_statement(&self, prepared_statement: &Statement) {
        prepared_statement.execute(&[&self.peptide_id, &self.protein_id]);
    }
}