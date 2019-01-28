use proteomic::models::enzyms::digest_enzym::DigestEnzym;
use proteomic::models::persistable::{Persistable, handle_postgres_error};
use proteomic::models::peptides::peptide::Peptide;
use proteomic::models::peptide_protein_association::PeptideProteinAssociation;


// attributes max_number_of_missed_cleavages should be unsigned, but it will passed to Peptides. see comment for Peptide for additional reasons
pub struct Trypsin<'t> {
    name: &'static str,
    shortcut: &'static str, // the returnes string should only 5 charachters long because field 'digest_enzym' of table peptides in databse is limited to 5 characters
    digist_regex: onig::Regex,
    digest_replace: &'static str,
    max_number_of_missed_cleavages: i16,
    // replace the next two attributes with a range in the feature: https://doc.rust-lang.org/std/ops/struct.Range.html#method.contains
    min_peptide_length: usize,
    max_peptide_length: usize,
    database_connection: &'t postgres::Connection,
    peptide_create_statement: postgres::stmt::Statement<'t>,
    peptide_exists_statement: postgres::stmt::Statement<'t>,
    pp_association_create_statement: postgres::stmt::Statement<'t>,
    pp_association_exists_statement: postgres::stmt::Statement<'t>,
}

impl<'t> DigestEnzym<'t> for Trypsin<'t> {
    fn new(database_connection: &'t postgres::Connection, max_number_of_missed_cleavages: i16, min_peptide_length: usize, max_peptide_length: usize) -> Self {
        Self {
            name: "Trypsin",
            shortcut: "try",
            digist_regex: onig::Regex::new(r"(?<=[KR])(?!P)").unwrap(),
            digest_replace: "$before $after",
            max_number_of_missed_cleavages: max_number_of_missed_cleavages,
            min_peptide_length: min_peptide_length,
            max_peptide_length: max_peptide_length,
            database_connection: database_connection,
            peptide_create_statement: match database_connection.prepare_cached(Peptide::create_query()) {
                Ok(statement) => statement,
                Err(err) => panic!("proteomic::models::enzyms::trypsin::Trypsin.create_database_connection_and_prepared_statements(): Error at @: {}", handle_postgres_error(&err))
            },
            peptide_exists_statement: match database_connection.prepare_cached(Peptide::exists_query()) {
                Ok(statement) => statement,
                Err(err) => panic!("proteomic::models::enzyms::trypsin::Trypsin.create_database_connection_and_prepared_statements(): Error at @: {}", handle_postgres_error(&err))
            },
            pp_association_create_statement: match database_connection.prepare_cached(PeptideProteinAssociation::create_query()) {
                Ok(statement) => statement,
                Err(err) => panic!("proteomic::models::enzyms::trypsin::Trypsin.create_database_connection_and_prepared_statements(): Error at @: {}", handle_postgres_error(&err))
            },
            pp_association_exists_statement: match database_connection.prepare_cached(PeptideProteinAssociation::exists_query()) {
                Ok(statement) => statement,
                Err(err) => panic!("proteomic::models::enzyms::trypsin::Trypsin.create_database_connection_and_prepared_statements(): Error at @: {}", handle_postgres_error(&err))
            }
        }
    }

    fn get_name(&self) -> &str {
        return self.name;
    }

    fn get_shortcut(&self) -> &str {
        return self.shortcut;
    }

    fn get_max_number_of_missed_cleavages(&self) -> i16 {
        return self.max_number_of_missed_cleavages;
    }

    fn get_digest_regex(&self) -> &onig::Regex {
        return &self.digist_regex;
    }

    fn get_digest_replace(&self) -> &'static str {
        return self.digest_replace;
    }

    fn get_min_peptide_length(&self) -> usize {
        return self.min_peptide_length;
    }

    fn get_max_peptide_length(&self) -> usize {
        return self.max_peptide_length;
    }

    fn get_database_connection(&self) -> &postgres::Connection {
        return self.database_connection;
    }

    fn get_peptide_create_statement(&self) -> &postgres::stmt::Statement {
        return &self.peptide_create_statement;
    }

    fn get_peptide_exists_statement(&self) -> &postgres::stmt::Statement {
        return &self.peptide_exists_statement;
    }

    fn get_pp_association_create_statement(&self) -> &postgres::stmt::Statement {
        return &self.pp_association_create_statement;
    }

    fn get_pp_association_exists_statement(&self) -> &postgres::stmt::Statement {
        return &self.pp_association_exists_statement;
    }
}