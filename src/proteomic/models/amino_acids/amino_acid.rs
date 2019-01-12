use proteomic::models::mass;
use proteomic::models::mass::neutral_loss::NeutralLoss;

const A: (&'static str, char, &'static str, &'static str, f64, f64, f32) = ("Alanine", 'A', "Ala", "C3H5ON", 71.03711, 71.0788, 9.1);
/** B is an average between N and D */
const B: (&'static str, char, &'static str, &'static str, f64, f64, f32) = ("Asparagine or aspartic acid", 'B', "Asx", "", 114.53495, 114.5962, 9.2);
const R: (&'static str, char, &'static str, &'static str, f64, f64, f32) = ("Arginine", 'R', "Arg", "C6H12ON4", 156.10111, 156.1875, 5.7);
const N: (&'static str, char, &'static str, &'static str, f64, f64, f32) = ("Asparagine", 'N', "Asn", "C4H6O2N2", 114.04293, 114.1038, 3.8);
const D: (&'static str, char, &'static str, &'static str, f64, f64, f32) = ("Aspartic acid", 'D', "Asp", "C4H5O3N", 115.02694, 115.0886, 5.4);
const C: (&'static str, char, &'static str, &'static str, f64, f64, f32) = ("Cysteine", 'C', "Cys", "C3H5ONS", 103.00919, 103.1388, 1.1);
const E: (&'static str, char, &'static str, &'static str, f64, f64, f32) = ("Glutamic acid", 'E', "Glu", "C5H7O3N", 129.04259, 129.1155, 6.1);
const Q: (&'static str, char, &'static str, &'static str, f64, f64, f32) = ("Glutamine", 'Q', "Gln", "C5H8O2N2", 128.05858, 128.1307, 3.7);
const G: (&'static str, char, &'static str, &'static str, f64, f64, f32) = ("Glycine", 'G', "Gly", "C2H3ON", 57.02146, 57.0519, 7.3);
const H: (&'static str, char, &'static str, &'static str, f64, f64, f32) = ("Histidine", 'H', "His", "C6H7ON3", 137.05891, 137.1411, 2.1);
const I: (&'static str, char, &'static str, &'static str, f64, f64, f32) = ("Isoleucine", 'I', "Ile", "C6H11ON", 113.08406, 113.1594, 5.6);
const L: (&'static str, char, &'static str, &'static str, f64, f64, f32) = ("Leucine", 'L', "Leu", "C6H11ON", 113.08406, 113.1594, 9.8);
const J: (&'static str, char, &'static str, &'static str, f64, f64, f32) = ("Isoleucine or Leucine", 'J', "Ile or Leu", "C6H11ON", 113.08406, 113.1594, 15.4);
const K: (&'static str, char, &'static str, &'static str, f64, f64, f32) = ("Lysine", 'K', "Lys", "C6H12ON2", 128.09496, 128.1741, 4.9);
const M: (&'static str, char, &'static str, &'static str, f64, f64, f32) = ("Methionine", 'M', "Met", "C5H9ONS", 131.04049, 131.1926, 2.3);
const F: (&'static str, char, &'static str, &'static str, f64, f64, f32) = ("Phenylalanine", 'F', "Phe", "C9H9ON", 147.06841, 147.1766, 3.9);
const P: (&'static str, char, &'static str, &'static str, f64, f64, f32) = ("Proline", 'P', "Pro", "C5H7ON", 97.05276, 97.1167, 4.8);
const O: (&'static str, char, &'static str, &'static str, f64, f64, f32) = ("Pyrrolysine", 'O', "Pyl", "C12H21N3O3", 109.0528, 109.13, 0.0);
const S: (&'static str, char, &'static str, &'static str, f64, f64, f32) = ("Serine", 'S', "Ser", "C3H5O2N", 87.03203, 87.0782, 6.6);
const T: (&'static str, char, &'static str, &'static str, f64, f64, f32) = ("Threonine", 'T', "Thr", "C4H7O2N", 101.04768, 101.1051, 5.5);
const U: (&'static str, char, &'static str, &'static str, f64, f64, f32) = ("Selenocysteine", 'U', "SeC", "C3H5NOSe", 150.95363, 150.0379, 0.0);
const V: (&'static str, char, &'static str, &'static str, f64, f64, f32) = ("Valine", 'V', "Val", "C5H9ON", 99.06841, 99.1326, 6.8);
const W: (&'static str, char, &'static str, &'static str, f64, f64, f32) = ("Tryptophan", 'W', "Trp", "C11H10ON2", 186.07931, 186.2132, 1.3);
/** Some Search Engines and Databases used the X Amino Acid for unknown amino acids*/
const X: (&'static str, char, &'static str, &'static str, f64, f64, f32) = ("Unknown Amino Acid", 'X', "Xaa", "Unknown", 0.0, 0.0, 0.0);
const Y: (&'static str, char, &'static str, &'static str, f64, f64, f32) = ("Tyrosine", 'Y', "Tyr", "C9H9O2N", 163.06333, 163.1760, 2.9);
/** Z is an average btween E and Q */
const Z: (&'static str, char, &'static str, &'static str, f64, f64, f32) = ("Glutamine or glutamic acid", 'Z', "Glx", "", 128.55059, 128.6231, 9.8);


pub struct AminoAcid {
    name: &'static str,
    one_letter_code: char,
    three_letter_code: &'static str,
    chemical_formula: &'static str,
    mono_mass: i64,
    average_mass: i64,
    distribution: i8
}

impl AminoAcid {
    pub fn new(amino_acid_tupel: (&'static str, char, &'static str, &'static str, f64, f64, f32)) -> AminoAcid {
        return AminoAcid {
            name: amino_acid_tupel.0,
            one_letter_code: amino_acid_tupel.1,
            three_letter_code: amino_acid_tupel.2,
            chemical_formula: amino_acid_tupel.3,
            mono_mass: mass::convert_mass_to_int(amino_acid_tupel.4),
            average_mass: mass::convert_mass_to_int(amino_acid_tupel.5),
            distribution: (amino_acid_tupel.6 * 10.0) as i8
        }
    }

    pub fn get_name (&self) -> &'static str {
        return self.name;
    }
    pub fn get_one_letter_code (&self) -> char {
        return self.one_letter_code;
    }
    pub fn get_three_letter_code (&self) -> &'static str {
        return self.three_letter_code;
    }
    pub fn get_chemical_formular (&self) -> &'static str {
        return self.chemical_formula;
    }
    pub fn get_mono_mass (&self) -> i64 {
        return self.mono_mass;
    }

    pub fn get_average_mass (&self) -> i64 {
        return self.average_mass;
    }

    pub fn get_distribution(&self) -> i8 {
        return self.distribution;
    }

    pub fn get(one_letter_code: char) -> AminoAcid {
        match one_letter_code {
            'A' => return AminoAcid::new(A),
            'B' => return AminoAcid::new(B),
            'R' => return AminoAcid::new(R),
            'N' => return AminoAcid::new(N),
            'D' => return AminoAcid::new(D),
            'C' => return AminoAcid::new(C),
            'E' => return AminoAcid::new(E),
            'Q' => return AminoAcid::new(Q),
            'G' => return AminoAcid::new(G),
            'H' => return AminoAcid::new(H),
            'I' => return AminoAcid::new(I),
            'L' => return AminoAcid::new(L),
            'J' => return AminoAcid::new(J),
            'K' => return AminoAcid::new(K),
            'M' => return AminoAcid::new(M),
            'F' => return AminoAcid::new(F),
            'P' => return AminoAcid::new(P),
            'O' => return AminoAcid::new(O),
            'S' => return AminoAcid::new(S),
            'T' => return AminoAcid::new(T),
            'U' => return AminoAcid::new(U),
            'V' => return AminoAcid::new(V),
            'W' => return AminoAcid::new(W),
            'X' => return AminoAcid::new(X),
            'Y' => return AminoAcid::new(Y),
            'Z' => return AminoAcid::new(Z),
            _ => return AminoAcid::new(X)
        }
    }

    pub fn get_haviest() -> AminoAcid {
        return AminoAcid::new(W);
    }

    pub fn get_lightest() -> AminoAcid {
        return AminoAcid::new(G);
    }

    pub fn get_sequence_weight(sequence: &str) -> i64 {
        let mut weight: i64 = NeutralLoss::get("H2O").get_mono_mass();
        for amino_acid_one_letter_code in sequence.chars() {
            weight += AminoAcid::get(amino_acid_one_letter_code).get_mono_mass();
        }
        return weight;
    }

    pub fn gerneralize_sequence(sequence: &str) -> String {
        return sequence.replace("I", "J").replace("L", "J");
    }
}



