use std::env;
use std::fs::File;
use std::io::BufReader;
use std::io::prelude::*;

mod bio;

use bio::enzym::{DigstEnzym, Trypsin};
use bio::protein::Protein;


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
        if !string_line.starts_with(">") {
            aa_sequence.push_str(&string_line);
        } else {
            if accession.len() > 0 {
                let mut protein: Protein = Protein::new(accession, aa_sequence);
                let mut trypsin: Trypsin = Trypsin::new(3, 1, 50);
                trypsin.digest(&mut protein);
                //aa_sequence = String::new();
                break;
            }
            accession = string_line;
        }
    }
}

// fn main (){
//     let peps = vec!["A","B","C","D","E","F","G","H"];
//     let msc = vec![2, 3, 7];
//     let mut new_pep = String::new();
//     for i in 0..msc.len() {
//         new_pep.push_str(peps[msc[i]]);
//         if (i + 1) < (msc.len()) && msc[i+1] != msc[i] + 1 {
//             // add to peptides
//             new_pep = String::new();
//         }
//     }
//     // add to peptides
// }
