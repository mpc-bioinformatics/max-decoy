extern crate postgres;

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
            peptide_id: *peptide.get_primary_key(),
            protein_id: *protein.get_primary_key(),
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
    fn get_primary_key(&self) -> &(i32, i32) {
        return &(self.peptide_id, self.protein_id);
    }

    fn find(conn: &Connection, primary_key: &(i32, i32)) -> Result<Self, &'static str> {
        match conn.query("SELECT * FROM peptides_proteins WHERE peptide_id = $1 AND protein_id = $2 LIMIT 1", &[&primary_key.0, &primary_key.1]) {
            Ok(rows) =>{
                if rows.len() > 0 {
                    Ok(
                        PeptideProteinAssociation {
                            peptide_id: rows.get(0).get(0),
                            protein_id: rows.get(0).get(1),
                            is_persisted: true
                        }
                    )
                } else {
                    Err("protein-peptide-association not found")
                }
            },
            Err(_err) => Err("some error occured: PeptideProteinAssociation::find")
        }
    }

    fn find_by_unique_identifier(conn: &Connection, unique_identifier: &(i32, i32)) -> Result<Self, &'static str> {
        return Self::find(conn, &unique_identifier);
    }


    fn create(&self, conn: &postgres::Connection) -> bool {
        match conn.query(
            Self::get_insert_query(),
            &[&self.peptide_id, &self.protein_id]
        ) {
            Ok(rows) => {
                if rows.len() > 0 {
                    self.peptide_id = rows.get(0).get(0);
                    self.protein_id = rows.get(0).get(1);
                    self.is_persisted = true;
                    return true;
                } else {
                    // zero rows means there are a conflict on update, so the peptides exists already
                    match Self::find_by_unique_identifier(conn, &(self.peptide_id, self.protein_id)) {
                        Ok(peptide_protein_association) => {
                            self.peptide_id = peptide_protein_association.get_primary_key().0;
                            self.protein_id = peptide_protein_association.get_primary_key().1;
                            self.is_persisted = true;
                            return true;
                        },
                        Err(_err) => {
                            println!("ERROR: something is really odd with peptide_protein_association {:?}:\n\tcan't create it nor find it\n", (self.peptide_id, self.protein_id));
                            return false;
                        }
                    }
                }
            }
            Err(_err) => return false
        }
    }

    fn update(&self, conn: &postgres::Connection) -> bool {
        match conn.query(
            Self::get_update_query(),
            &[&self.peptide_id, &self.protein_id]
        ) {
            Ok(rows) => return rows.len() > 0,
            Err(_err) => return false
        }
    }

    fn save(&self, conn: &postgres::Connection) -> bool {
        if !self.is_persisted() {
            return self.update(conn);
        } else {
            return self.create(conn);
        }
    }


    fn get_insert_query() -> &'static str {
        return "INSERT INTO peptides_proteins (peptide_id, protein_id) VALUES ($1, $2) ON CONFLICT (peptide_id, protein_id) DO NOTHING RETURNING *";
    }

    fn get_update_query() -> &'static str {
        return "UPDATE peptides_proteins SET peptide_id = $1 AND protein_id = $2 WHERE peptide_id = $1 AND protein_id = $2";
    }


    fn exec_insert_statement(&self, prepared_statement: &postgres::stmt::Statement) -> bool {
        match prepared_statement.query(&[&self.peptide_id, &self.protein_id]) {
            Ok(rows) => {
                if rows.len() > 0 {
                    self.peptide_id = rows.get(0).get(0);
                    self.protein_id = rows.get(0).get(1);
                    self.is_persisted = true;
                    return true;
                } else {
                    return false
                }
            },
            Err(_err) => return false
        }
    }

    fn exec_update_statement(&self, prepared_statement: &postgres::stmt::Statement) -> bool {
        match prepared_statement.query(&[&self.peptide_id, &self.protein_id]) {
            Ok(rows) => return rows.len() > 0,
            Err(_err) => return false
        }
    }


    fn is_persisted(&self) -> bool {
        return self.is_persisted;
    }

}