pub struct AminoAcid {
    name: &'static str,
    one_letter_code: char,
    three_letter_code: &'static str,
    chemical_formula: &'static str,
    mono_mass: f64,
    average_mass: f64
}

impl AminoAcid {
    pub fn new(amino_acid_tupel: (&'static str, char, &'static str, &'static str, f64, f64)) -> AminoAcid {
        return AminoAcid {
            name: amino_acid_tupel.0,
            one_letter_code: amino_acid_tupel.1,
            three_letter_code: amino_acid_tupel.2,
            chemical_formula: amino_acid_tupel.3,
            mono_mass: amino_acid_tupel.4,
            average_mass: amino_acid_tupel.5  
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
    pub fn get_mono_mass (&self) -> f64 {
        return self.mono_mass;
    }
    pub fn get_average_mass (&self) -> f64 {
        return self.average_mass;
    }
}

const A: (&'static str, char, &'static str, &'static str, f64, f64) = ("Alanine", 'A', "Ala", "C3H5ON", 71.03711, 71.0788);
/** B is an average between N and D */
const B: (&'static str, char, &'static str, &'static str, f64, f64) = ("Asparagine or aspartic acid", 'B', "Asx", "", 114.53495, 114.5962);
const R: (&'static str, char, &'static str, &'static str, f64, f64) = ("Arginine", 'R', "Arg", "C6H12ON4", 156.10111, 156.1875);
const N: (&'static str, char, &'static str, &'static str, f64, f64) = ("Asparagine", 'N', "Asn", "C4H6O2N2", 114.04293, 114.1038);
const D: (&'static str, char, &'static str, &'static str, f64, f64) = ("Aspartic acid", 'D', "Asp", "C4H5O3N", 115.02694, 115.0886);
const C: (&'static str, char, &'static str, &'static str, f64, f64) = ("Cysteine", 'C', "Cys", "C3H5ONS", 103.00919, 103.1388);
const E: (&'static str, char, &'static str, &'static str, f64, f64) = ("Glutamic acid", 'E', "Glu", "C5H7O3N", 129.04259, 129.1155);
const Q: (&'static str, char, &'static str, &'static str, f64, f64) = ("Glutamine", 'Q', "Gln", "C5H8O2N2", 128.05858, 128.1307);
const G: (&'static str, char, &'static str, &'static str, f64, f64) = ("Glycine", 'G', "Gly", "C2H3ON", 57.02146, 57.0519);
const H: (&'static str, char, &'static str, &'static str, f64, f64) = ("Histidine", 'H', "His", "C6H7ON3", 137.05891, 137.1411);
const I: (&'static str, char, &'static str, &'static str, f64, f64) = ("Isoleucine", 'I', "Ile", "C6H11ON", 113.08406, 113.1594);
const L: (&'static str, char, &'static str, &'static str, f64, f64) = ("Leucine", 'L', "Leu", "C6H11ON", 113.08406, 113.1594);
const J: (&'static str, char, &'static str, &'static str, f64, f64) = ("Isoleucine or Leucine", 'J', "Ile or Leu", "C6H11ON", 113.08406, 113.1594);
const K: (&'static str, char, &'static str, &'static str, f64, f64) = ("Lysine", 'K', "Lys", "C6H12ON2", 128.09496, 128.1741);
const M: (&'static str, char, &'static str, &'static str, f64, f64) = ("Methionine", 'M', "Met", "C5H9ONS", 131.04049, 131.1926);
const F: (&'static str, char, &'static str, &'static str, f64, f64) = ("Phenylalanine", 'F', "Phe", "C9H9ON", 147.06841, 147.1766);
const P: (&'static str, char, &'static str, &'static str, f64, f64) = ("Proline", 'P', "Pro", "C5H7ON", 97.05276, 97.1167);
const O: (&'static str, char, &'static str, &'static str, f64, f64) = ("Pyrrolysine", 'O', "Pyl", "C12H21N3O3", 109.0528, 109.13);
const S: (&'static str, char, &'static str, &'static str, f64, f64) = ("Serine", 'S', "Ser", "C3H5O2N", 87.03203, 87.0782);
const T: (&'static str, char, &'static str, &'static str, f64, f64) = ("Threonine", 'T', "Thr", "C4H7O2N", 101.04768, 101.1051);
const U: (&'static str, char, &'static str, &'static str, f64, f64) = ("Selenocysteine", 'U', "SeC", "C3H5NOSe", 150.95363, 150.0379);
const V: (&'static str, char, &'static str, &'static str, f64, f64) = ("Valine", 'V', "Val", "C5H9ON", 99.06841, 99.1326);
const W: (&'static str, char, &'static str, &'static str, f64, f64) = ("Tryptophan", 'W', "Trp", "C11H10ON2", 186.07931, 186.2132);
/** Some Search Engines and Daabases used the X Amino Acid for unknown amino acids*/
const X: (&'static str, char, &'static str, &'static str, f64, f64) = ("Unknown Amino Acid", 'X', "Xaa", "Unknown", 0.0, 0.0);
const Y: (&'static str, char, &'static str, &'static str, f64, f64) = ("Tyrosine", 'Y', "Tyr", "C9H9O2N", 163.06333, 163.1760);
/** Z is an average btween E and Q */
const Z: (&'static str, char, &'static str, &'static str, f64, f64) = ("Glutamine or glutamic acid", 'Z', "Glx", "", 128.55059, 128.6231);


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