use std::env;
use std::fs::File;
use std::io::BufReader;
use std::io::prelude::*;

mod proteomic;

use proteomic::utility::enzym::{DigestEnzym, Trypsin};
use proteomic::models::protein::Protein;


fn main() {
    let args: Vec<String> = env::args().collect();

    let filename = &args[1];

    // --snip--
    println!("In file {}", filename);

    let fasta_file = File::open(filename).expect("fasta file not found");
    let fasta_file = BufReader::new(fasta_file);

    let mut accession: String = String::new();
    let mut aa_sequence = String::new();

    for line in fasta_file.lines() {
        let string_line = line.unwrap();
        // strip and replace/error whitespaces
        if !string_line.starts_with(">") {
            aa_sequence.push_str(&string_line);
        } else {
            if accession.len() > 0 {
                let mut protein: Protein = Protein::new(accession, aa_sequence);
                let mut trypsin: Trypsin = Trypsin::new(3, 1, 50);
                let peptides = trypsin.digest(&mut protein);
                for pep in peptides {
                    println!("{}", pep);
                }
                //aa_sequence = String::new();
                break;
            }
            accession = string_line;
        }
    }
    // add last line
}
