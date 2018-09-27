extern crate postgres;

use std::hash::{Hash, Hasher};

use self::postgres::Connection;
use self::postgres::rows::Rows;
use self::postgres::Result;
use self::postgres::stmt::Statement;

use proteomic::models::peptide::Peptide;
use proteomic::models::protein::Protein;
use proteomic::models::persistable::Persistable;
use proteomic::models::collection::Collectable;

pub struct PeptideProteinAssociation {
    // protein_id: i32,
    // peptide_id: i32,
    peptide_aa_sequence: String,
    protein_accession: String,
    collection_identifier: String
}

impl PeptideProteinAssociation {
    pub fn new(peptide: &Peptide, protein: &Protein) -> PeptideProteinAssociation {
        return PeptideProteinAssociation {
            // protein_id: -1,
            // peptide_id: -1,
            peptide_aa_sequence: peptide.get_aa_sequence().clone(),
            protein_accession: protein.get_accession().clone(),
            collection_identifier: PeptideProteinAssociation::build_collection_identifier(peptide.get_aa_sequence(), protein.get_accession())
        }
    }

    pub fn print(&self) {
        println!("{} = {} <=> {}",self.collection_identifier, self.peptide_aa_sequence, self.protein_accession);
    }

    fn build_collection_identifier(peptide_aa_sequence: &str, protein_accession: &str) -> String {
        let mut identifier:String = peptide_aa_sequence.to_owned();
        identifier.push_str(protein_accession);
        return identifier;
    }

    // pub fn get_actual_ids(&mut self) {
    //     let conn = super::get_db_connection();
    //     for row in &conn.query("SELECT id FROM peptides WHERE aa_sequence = $1", &[&self.peptide_aa_sequence]).unwrap() {
    //         self.peptide_id = row.get("id");
    //         break;
    //     }
    //     for row in &conn.query("SELECT id FROM proteins WHERE accession = $1", &[&self.protein_accession]).unwrap() {
    //         self.protein_id = row.get("id");
    //         break;
    //     }
    // }

    pub fn get_peptide_aa_sequence(&self) -> &String {
        return &self.peptide_aa_sequence;
    }

    pub fn get_protein_accession(&self) -> &String {
        return &self.protein_accession;
    }

    pub fn create(&mut self, conn: &Connection) {
        self.execute_insert_query(&conn);
    }


    pub fn exists(&self, conn: &Connection) -> bool {
        for row in self.exists_query(&conn).unwrap().iter() {
            return row.get::<usize, bool>(0);
        }
        return false;
    }
}

impl Persistable for PeptideProteinAssociation {
    fn get_id(&self) -> i32 {
        return 0;
    }

    fn get_insert_statement() -> &'static str {
        return "INSERT INTO peptides_proteins (peptide_aa_sequence, protein_accession) VALUES ($1, $2) ON CONFLICT DO NOTHING";
        // return "INSERT INTO peptides_proteins (peptide_id, protein_id) VALUES ($1, $2) ON CONFLICT DO NOTHING RETURNING id";
    }

    fn execute_insert_query(&self, connection: &Connection) -> Result<Rows> {
        return connection.query(
            PeptideProteinAssociation::get_insert_statement(),
            &[&self.peptide_aa_sequence, &self.protein_accession]
            // &[&self.peptide_id, &self.protein_id]
        );
    }

    fn execute_insert_statement(&self, prepared_statement: &Statement) {
        prepared_statement.execute(&[&self.peptide_aa_sequence, &self.protein_accession]);
        // prepared_statement.execute(&[&self.peptide_id, &self.protein_id]);
    }

    // update does not make sens yet, see 'where'-clause
    fn get_update_statement() -> &'static str {
        return "UPDATE peptides_proteins SET peptide_aa_sequence = $1, protein_id = $2 WHERE protein_id = 0 AND protein_id = 0";
        // return "UPDATE peptides_proteins SET peptide_id = $1, protein_id = $2 WHERE protein_id = 0 AND protein_id = 0";
    }

    fn execute_update_query(&self, connection: &Connection) -> Result<Rows> {
        return connection.query(
            PeptideProteinAssociation::get_update_statement(),
            &[&self.peptide_aa_sequence, &self.protein_accession]
            // &[&self.peptide_id, &self.protein_id]
        );
    }

    fn execute_update_statement(&self, prepared_statement: &Statement) {
        prepared_statement.execute(&[&self.peptide_aa_sequence, &self.protein_accession]);
        // prepared_statement.execute(&[&self.peptide_id, &self.protein_id]);
    }

    fn get_exists_statement() -> &'static str {
        return "SELECT EXISTS(SELECT 1 FROM peptides_proteins WHERE peptide_id = $1 AND protein_id = $2)";
    }

    fn exists_query(&self, connection: &Connection) -> Result<Rows> {
        return connection.query(
            PeptideProteinAssociation::get_exists_statement(),
            &[&self.peptide_aa_sequence, &self.protein_accession]
            // &[&self.peptide_id, &self.protein_id]
        );
    }
}

impl Collectable for PeptideProteinAssociation {
    fn get_collection_identifier(&self) -> &String {
      return &self.collection_identifier;
    }
}

// PartialEq-implementation to use this type in a HashSet
impl PartialEq for PeptideProteinAssociation {
    fn eq(&self, other: &PeptideProteinAssociation) -> bool {
       return (self.peptide_aa_sequence == *other.get_peptide_aa_sequence()) & (self.protein_accession == *other.get_protein_accession());
    }
}

// Eq-implementation to use this type in a HashSet
impl Eq for PeptideProteinAssociation {}

// Hash-implementation to use this type in a HashSet
impl Hash for PeptideProteinAssociation {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.collection_identifier.hash(state);
    }
}