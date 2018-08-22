use std::env;
use std::fs::File;
use std::io::BufReader;
use std::io::prelude::*;
use std::collections::HashSet;

mod proteomic;

use proteomic::utility::enzym::{DigestEnzym, Trypsin};
use proteomic::models::protein::Protein;
use proteomic::models::peptide::Peptide;

fn main() {
    let args: Vec<String> = env::args().collect();

    let filename = &args[1];

    // --snip--
    println!("In file {}", filename);

    let fasta_file = File::open(filename).expect("fasta file not found");
    let fasta_file = BufReader::new(fasta_file);

    let mut header: String = String::new();
    let mut aa_sequence = String::new();
    let trypsin: Trypsin = Trypsin::new(3, 1, 50);

    let mut peptides: HashSet<Peptide> = HashSet::new();

    for line in fasta_file.lines() {
        // trim
        let string_line = line.unwrap().as_mut_str().trim().to_owned()  ;
        if !string_line.starts_with(">") {
            aa_sequence.push_str(&string_line);
        } else {
            if header.len() > 0 {
                let mut protein: Protein = Protein::new(header, aa_sequence);
                // throw error if aa sequence has whitespaces
                peptides.extend(trypsin.digest(&mut protein));
                // for pep in &peptides {
                //     pep.print();
                // }
                aa_sequence = String::new();
            }
            header = string_line;
        }
    }
    // process last protein
    let mut protein: Protein = Protein::new(header, aa_sequence);
    // throw error if aa sequence has whitespaces
    peptides.extend(trypsin.digest(&mut protein));
    println!("Peptides created: {}", peptides.len());
}
