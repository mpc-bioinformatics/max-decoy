pub mod digest_enzym;
pub mod trypsin;
pub mod digest_summary;
pub mod transaction_summary;

use self::digest_enzym::DigestEnzym;
use self::trypsin::Trypsin;

pub fn get<'t>(enzym_name: &str, database_connection: &'t postgres::Connection, max_number_of_missed_cleavages: u8, min_peptide_length: usize, max_peptide_length: usize) -> impl DigestEnzym<'t> {
    match enzym_name.to_lowercase().as_str() {
        "trypsin" => Trypsin::new(database_connection, max_number_of_missed_cleavages, min_peptide_length, max_peptide_length),
        _ => Trypsin::new(database_connection, max_number_of_missed_cleavages, min_peptide_length, max_peptide_length)
    }
}