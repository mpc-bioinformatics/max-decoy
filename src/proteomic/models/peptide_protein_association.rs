extern crate postgres;

use std::error::Error;

use self::postgres::Connection;

use proteomic::models::peptide::Peptide;
use proteomic::models::protein::Protein;
use proteomic::models::persistable::Persistable;

pub struct PeptideProteinAssociation {
    peptide_id: i32,
    protein_id: i32,
    is_persisted: bool
}

impl PeptideProteinAssociation {
    pub fn new(peptide: &Peptide, protein: &Protein) -> PeptideProteinAssociation {
        return PeptideProteinAssociation {
            peptide_id: peptide.get_primary_key(),
            protein_id: protein.get_primary_key(),
            is_persisted: false

        }
    }

    pub fn get_peptide_id(&self) -> i32 {
        return self.peptide_id;
    }

    pub fn get_protein_id(&self) -> i32 {
        return self.protein_id;
    }
}

impl Persistable<PeptideProteinAssociation, (i32, i32), (i32, i32)> for PeptideProteinAssociation {
    fn get_primary_key(&self) -> (i32, i32) {
        return (self.peptide_id, self.protein_id);
    }

    fn find(conn: &Connection, primary_key: &(i32, i32)) -> Result<Self, String> {
        return Self::find_by_unique_identifier(conn, primary_key);
    }

    fn find_by_unique_identifier(conn: &Connection, unique_identifier: &(i32, i32)) -> Result<Self, String> {
        match conn.query(
            "SELECT * FROM peptides WHERE peptide_id = $1 and protein_id = $2 LIMIT 1",
            &[&unique_identifier.0, &unique_identifier.1]
        ) {
            Ok(ref rows) if rows.len() > 0 =>  Ok(
                PeptideProteinAssociation {
                    peptide_id: rows.get(0).get(0),
                    protein_id: rows.get(0).get(1),
                    is_persisted: true
                }
            ),
            Ok(_rows) => Err("protein-peptide-association not found".to_owned()),
            Err(err) => Err(err.description().to_owned())
        }
    }


    fn create(&mut self, conn: &postgres::Connection) -> Result<(), String> {
        match conn.query(
            Self::get_insert_query(),
            &[&self.peptide_id, &self.protein_id]
        ) {
            Ok(ref rows) if rows.len() > 0 => {
                self.peptide_id = rows.get(0).get(0);
                self.protein_id = rows.get(0).get(1);
                self.is_persisted = true;
                return Ok(());
            }
            Ok(_rows) => {
                // zero rows means there are a conflict on update, so the peptides exists already
                match Self::find_by_unique_identifier(conn, &(self.peptide_id, self.protein_id)) {
                    Ok(peptide_protein_association) => {
                        self.peptide_id = peptide_protein_association.get_primary_key().0;
                        self.protein_id = peptide_protein_association.get_primary_key().1;
                        self.is_persisted = true;
                        return Ok(());
                    },
                    Err(_err) => Err(format!("cannot inser not find peptides-protein-accession '({}, {})'", self.peptide_id, self.protein_id))
                }
            }
            Err(err) => Err(err.description().to_owned())
        }
    }

    fn update(&mut self, conn: &postgres::Connection) -> Result<(), String> {
        match conn.query(
            Self::get_update_query(),
            &[&self.peptide_id, &self.protein_id]
        ) {
            Ok(ref rows) if rows.len() > 0 => Ok(()),
            Ok(_rows) => Err("updateing peptide-protein-association does not return anything".to_owned()),
            Err(err) => Err(err.description().to_owned())
        }
    }

    fn save(&mut self, conn: &postgres::Connection) -> Result<(), String> {
        if self.is_persisted() {
            return self.update(conn);
        } else {
            return self.create(conn);
        }
    }


    fn get_select_primary_key_by_unique_identifier_query() -> &'static str {
        return "SELECT peptide_id, protein_id FROM peptides_proteins WHERE peptide_id = $1 and protein_id = $2 LIMIT 1";
    }

    fn get_insert_query() -> &'static str {
        return "INSERT INTO peptides_proteins (peptide_id, protein_id) VALUES ($1, $2) ON CONFLICT (peptide_id, protein_id) DO NOTHING RETURNING *";
    }

    fn get_update_query() -> &'static str {
        return "UPDATE peptides_proteins SET peptide_id = $1 AND protein_id = $2 WHERE peptide_id = $1 AND protein_id = $2";
    }


    fn exec_select_primary_key_by_unique_identifier_statement(&mut self, prepared_statement: &postgres::stmt::Statement) -> Result<(), String> {
        match prepared_statement.query(&[&self.peptide_id, &self.protein_id]) {
            Ok(ref rows) if rows.len() > 0 => {
                self.peptide_id = rows.get(0).get(0);
                self.protein_id = rows.get(0).get(1);
                self.is_persisted = true;
                return Ok(());
            },
            Ok(_rows) => Err("peptide-protein-association not found".to_owned()),
            Err(err) => Err(err.description().to_owned())
        }
    }

    fn exec_insert_statement(&mut self, prepared_statement: &postgres::stmt::Statement) -> Result<(), String> {
        match prepared_statement.query(&[&self.peptide_id, &self.protein_id]) {
            Ok(ref rows) if rows.len() > 0 => {
                self.peptide_id = rows.get(0).get(0);
                self.protein_id = rows.get(0).get(1);
                self.is_persisted = true;
                return Ok(());
            },
            Ok(_rows) => Err("inseting peptide-protein-association does not return anything".to_owned()),
            Err(err) => Err(err.description().to_owned())
        }
    }

    fn exec_update_statement(&mut self, prepared_statement: &postgres::stmt::Statement) -> Result<(), String> {
        match prepared_statement.query(&[&self.peptide_id, &self.protein_id]) {
            Ok(ref rows) if rows.len() > 0 => Ok(()),
            Ok(_rows) => Err("updating peptide-protein-association does not return anything".to_owned()),
            Err(err) => Err(err.description().to_owned())
        }
    }


    fn is_persisted(&self) -> bool {
        return self.is_persisted;
    }

}