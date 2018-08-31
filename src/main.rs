extern crate time;

use std::env;
use std::fs::File;
use std::io::BufReader;
use std::io::prelude::*;

mod proteomic;
use proteomic::utility::enzym::{DigestEnzym, Trypsin};
use proteomic::models::collection::Collection;
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
    let trypsin: Trypsin = Trypsin::new(2, 6, 50);

    let mut peptides: Collection<Peptide> = Collection::new();
    let mut proteins: Collection<Protein> = Collection::new();
    let mut overall_protein_counter: usize = 0;
    let mut overall_peptide_counter: usize = 0;

    let start_time: f64 = time::precise_time_s();
    for line in fasta_file.lines() {
        // trim
        let string_line = line.unwrap().as_mut_str().trim().to_owned()  ;
        if !string_line.starts_with(">") {
            aa_sequence.push_str(&string_line);
        } else {
            if header.len() > 0 {
                let mut protein: Protein = Protein::new(header.clone(), aa_sequence);
                // throw error if aa sequence has whitespaces
                trypsin.digest(&mut protein, &mut peptides);
                proteins.add(protein);
                aa_sequence = String::new();
                if peptides.len() == 100000 {
                    overall_protein_counter += proteins.len();
                    overall_peptide_counter += peptides.len();
                    proteins.save();
                    peptides.save();
                    proteins.clear();
                    peptides.clear();
                }
            }
            header = string_line;
        }
    }
    // process last protein
    let mut protein: Protein = Protein::new(header, aa_sequence);
    // throw error if aa sequence has whitespaces
    trypsin.digest(&mut protein, &mut peptides);
    overall_protein_counter += proteins.len();
    overall_peptide_counter += peptides.len();
    peptides.save();
    proteins.save();
    let stop_time: f64 = time::precise_time_s();
    println!("Proteins processed: {}", overall_protein_counter);
    println!("  Peptides created: {}", overall_peptide_counter);
    println!("Need {} s", (stop_time - start_time));
}