extern crate postgres;

use self::postgres::Connection;

use proteomic::models::persistable::{handle_postgres_error, Persistable, QueryError, QueryOk, FromSqlRowError};
use proteomic::models::peptide::Peptide;
use proteomic::models::protein::Protein;

pub struct PeptideProteinAssociation {
    peptide_id: i64,
    protein_id: i64
}

impl PeptideProteinAssociation {
    pub fn new(peptide: &Peptide, protein: &Protein) -> PeptideProteinAssociation {
        return PeptideProteinAssociation {
            peptide_id: peptide.get_primary_key(),
            protein_id: protein.get_primary_key()
        }
    }

    pub fn get_peptide_id(&self) -> i64 {
        return self.peptide_id;
    }

    pub fn get_protein_id(&self) -> i64 {
        return self.protein_id;
    }
}

impl Persistable<PeptideProteinAssociation, (i64, i64), (i64, i64)> for PeptideProteinAssociation {
    fn from_sql_row(row: &postgres::rows::Row) -> Result<Self, FromSqlRowError> {
        return Ok(
            Self {
                peptide_id: row.get(0),
                protein_id: row.get(1)
            }
        )
    }

    fn set_primary_key_from_sql_row(&mut self, row: &postgres::rows::Row) {
        self.peptide_id = row.get(0);
        self.protein_id = row.get(1);
    }

    fn invalidate_primary_key(&mut self) {
        self.peptide_id = 0;
        self.protein_id = 0;
    }

    fn get_primary_key(&self) -> (i64, i64) {
        return (self.peptide_id, self.protein_id);
    }

    fn get_table_name() -> &'static str {
        return "peptides_proteins";
    }

    fn is_persisted(&self) -> bool {
        return (self.peptide_id > 0) & (self.protein_id > 0);
    }

    fn find_query() -> &'static str {
        return "SELECT * FROM peptides_proteins WHERE peptide_id = $1 and protein_id = $2 LIMIT 1;";
    }

    fn create_query() -> &'static str {
        return "INSERT INTO peptides_proteins (peptide_id, protein_id) VALUES ($1, $2) ON CONFLICT (peptide_id, protein_id) DO NOTHING RETURNING *;";
    }

    fn create_attributes(&self) -> Box<Vec<&postgres::types::ToSql>>{
        return Box::new(vec![&self.peptide_id, &self.protein_id]);
    }

    fn update_query() -> &'static str{
        return "UPDATE peptides_proteins SET peptide_id = $1 AND protein_id = $2 WHERE peptide_id = $1 AND protein_id = $2;";
    }

    fn update_attributes(&self) -> Box<Vec<&postgres::types::ToSql>>{
        return Box::new(vec![&self.peptide_id, &self.protein_id]);
    }

    fn delete_query() -> &'static str {
        return "DELETE FROM peptides_proteins WHERE peptide_id = $1 AND protein_id = $2;";
    }

    fn delete_attributes(&self) -> Box<Vec<&postgres::types::ToSql>> {
        return Box::new(vec![&self.peptide_id, &self.protein_id]);
    }

    fn delete_all_query() -> &'static str {
        return "DELETE FROM peptides_proteins WHERE peptide_id IS NOT NULL AND protein_id IS NOT NULL;";
    }

    fn exists_query() -> &'static str {
        "SELECT * FROM peptides_proteins WHERE peptide_id = $1 AND protein_id = $2 LIMIT 1;"
    }

    fn exists_attributes(&self) -> Box<Vec<&postgres::types::ToSql>> {
        return Box::new(vec![&self.peptide_id, &self.protein_id]);
    }

    fn before_delete_hook(&self) -> Result<(), QueryError> {return Ok(());}
}
