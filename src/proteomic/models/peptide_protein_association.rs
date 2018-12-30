extern crate postgres;

use self::postgres::Connection;

use proteomic::models::peptide::Peptide;
use proteomic::models::protein::Protein;
use proteomic::models::persistable::Persistable;

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
    fn from_sql_row(row: &postgres::rows::Row) -> Result<Self, String> {
        return Ok(
            Self {
                peptide_id: row.get(0),
                protein_id: row.get(1)
            }
        )
    }

    fn get_primary_key(&self) -> (i64, i64) {
        return (self.peptide_id, self.protein_id);
    }

    fn find(conn: &Connection, primary_key: &(i64, i64)) -> Result<Self, String> {
        return Self::find_by_unique_identifier(conn, primary_key);
    }

    fn find_by_unique_identifier(conn: &Connection, unique_identifier: &(i64, i64)) -> Result<Self, String> {
        match conn.query(
            "SELECT * FROM peptides WHERE peptide_id = $1 and protein_id = $2 LIMIT 1",
            &[&unique_identifier.0, &unique_identifier.1]
        ) {
            Ok(ref rows) if rows.len() > 0 => Self::from_sql_row(&rows.get(0)),
            Ok(_rows) => Err("NOHIT".to_owned()),
            Err(err) => Err(err.code().unwrap().code().to_owned())
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
                return Ok(());
            }
            Ok(_rows) => {
                // zero rows means there are a conflict on update, so the peptides exists already
                match Self::find_by_unique_identifier(conn, &(self.peptide_id, self.protein_id)) {
                    Ok(peptide_protein_association) => {
                        self.peptide_id = peptide_protein_association.get_primary_key().0;
                        self.protein_id = peptide_protein_association.get_primary_key().1;
                        return Ok(());
                    },
                    Err(_err) => Err(format!("cannot inser not find peptides-protein-accession '({}, {})'", self.peptide_id, self.protein_id))
                }
            }
            Err(err) => Err(err.code().unwrap().code().to_owned())
        }
    }

    fn update(&mut self, conn: &postgres::Connection) -> Result<(), String> {
        match conn.query(
            Self::get_update_query(),
            &[&self.peptide_id, &self.protein_id]
        ) {
            Ok(ref rows) if rows.len() > 0 => Ok(()),
            Ok(_rows) => Err("NORET".to_owned()),
            Err(err) => Err(err.code().unwrap().code().to_owned())
        }
    }

    fn save(&mut self, conn: &postgres::Connection) -> Result<(), String> {
        if self.is_persisted() {
            return self.update(conn);
        } else {
            return self.create(conn);
        }
    }

    fn delete(&mut self, conn: &postgres::Connection) -> Result<(), String> {
        if !self.is_persisted() {
            return Err("ModifiedDecoys is not persisted".to_owned());
        }
        match conn.execute("DELETE FROM peptides_proteins WHERE peptide_id = $1 AND protein_id = $2;", &[&self.peptide_id, &self.protein_id]) {
            Ok(_) => {
                self.peptide_id = 0;
                self.protein_id = 0;
                return Ok(());
            },
            Err(err) => Err(format!("could not delete ModifiedDecoy from database; postgresql error is: {}", err))
        }
    }

    fn delete_all(conn: &postgres::Connection) -> Result<(), String> {
        match conn.execute("DELETE FROM peptides_proteins WHERE peptide_id IS NOT NULL AND protein_id IS NOT NULL;", &[]) {
            Ok(_) => Ok(()),
            Err(err) => Err(format!("could not delete ModifiedDecoys from database; postgresql error is: {}", err))
        }
    }

    fn select_where(conn: &postgres::Connection, conditions: &str, values: &[&postgres::types::ToSql]) -> Result<Vec<Self>, String> {
        let where_statement: String = format!("SELECT * FROM base_decoys WHERE {};", conditions);
        match conn.query(where_statement.as_str(), values) {
            Ok(ref rows) => {
                let records: Vec<Self> = Vec::new();
                for row in rows {
                    match Self::from_sql_row(&row) {
                        Ok(record) => records.push(record),
                        Err(err) => return Err(err)
                    }
                    
                }
                return Ok(records);
            },
            Err(err) => Err(format!("could not gether BaseDecoys from database; postgresql error is: {}", err))
        }
    }


    fn get_table_name() -> &'static str {
        return "peptides_proteins";
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
                return Ok(());
            },
            Ok(_rows) => Err("NORET".to_owned()),
            Err(err) => Err(err.code().unwrap().code().to_owned())
        }
    }

    fn exec_insert_statement(&mut self, prepared_statement: &postgres::stmt::Statement) -> Result<(), String> {
        match prepared_statement.query(&[&self.peptide_id, &self.protein_id]) {
            Ok(ref rows) if rows.len() > 0 => {
                self.peptide_id = rows.get(0).get(0);
                self.protein_id = rows.get(0).get(1);
                return Ok(());
            },
            Ok(_rows) => Err("NORET".to_owned()),
            Err(err) => Err(err.code().unwrap().code().to_owned())
        }
    }

    fn exec_update_statement(&mut self, prepared_statement: &postgres::stmt::Statement) -> Result<(), String> {
        match prepared_statement.query(&[&self.peptide_id, &self.protein_id]) {
            Ok(ref rows) if rows.len() > 0 => Ok(()),
            Ok(_rows) => Err("NORET".to_owned()),
            Err(err) => Err(err.code().unwrap().code().to_owned())
        }
    }


    fn is_persisted(&self) -> bool {
        return (self.peptide_id > 0) & (self.protein_id > 0);
    }
}
