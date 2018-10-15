pub mod digest_enzym;
pub mod trypsin;

use self::digest_enzym::DigestEnzym;
use self::trypsin::Trypsin;

pub fn get_enzym_by_name(enzym_name: &String, max_number_of_missed_cleavages: i16, min_peptide_length: usize, max_peptide_length: usize) -> impl DigestEnzym {
    match enzym_name.to_lowercase().as_str() {
        "trypsin" => Trypsin::new(max_number_of_missed_cleavages, min_peptide_length, max_peptide_length),
        _ => Trypsin::new(max_number_of_missed_cleavages, min_peptide_length, max_peptide_length)
    }
}