extern crate onig;

use proteomic::utility::enzyms::digest_enzym::DigestEnzym;

// attributes max_number_of_missed_cleavages should be unsigned, but it will passed to Peptides. see comment for Peptide for additional reasons
pub struct Trypsin {
    name: String,
    shortcut: String, // the returnes string should only 5 charachters long because field 'digest_enzym' of table peptides in databse is limited to 5 characters
    digist_regex: onig::Regex,
    digest_replace: &'static str,
    max_number_of_missed_cleavages: i16,
    // replace the next two attributes with a range in the feature: https://doc.rust-lang.org/std/ops/struct.Range.html#method.contains
    min_peptide_length: usize,
    max_peptide_length: usize
}

impl DigestEnzym for Trypsin {
    fn new (max_number_of_missed_cleavages: i16, min_peptide_length: usize, max_peptide_length: usize) -> Trypsin {
        Trypsin {
            name: "Trypsin".to_owned(),
            shortcut: "try".to_owned(),
            digist_regex: onig::Regex::new(r"(?<=[KR])(?!P)").unwrap(),
            digest_replace: "$before $after",
            max_number_of_missed_cleavages: max_number_of_missed_cleavages,
            min_peptide_length: min_peptide_length,
            max_peptide_length: max_peptide_length
        }
    }

    fn get_name(&self) -> &str {
        return self.name.as_str();
    }

    fn get_shortcut(&self) -> &str {
        return &self.shortcut;
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
}

impl Clone for Trypsin {
    fn clone(&self) -> Trypsin {
        return Trypsin::new(self.max_number_of_missed_cleavages, self.min_peptide_length, self.max_peptide_length);
    }
}

unsafe impl Send for Trypsin {}