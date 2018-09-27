extern crate time;

use std::env;
use std::fs::File;
use std::io::BufReader;
use std::io::prelude::*;

use std::collections::HashSet;

use proteomic::utility::enzym::{DigestEnzym, Trypsin};
use proteomic::models::collection::Collection;
use proteomic::models::protein::Protein;
use proteomic::models::peptide::Peptide;

#[test]
#[ignore]
fn test_digestion() {
    let mut path = env::current_dir().unwrap();
    path.push("fasta_files");
    path.push("tests");
    path.push("2018_05_20_UP000000625_escherichia_coli_reviewed");
    path.set_extension("fasta");

    // digest fasta file
    let fasta_file = File::open(path.to_str().unwrap()).expect("fasta file not found");
    let fasta_file = BufReader::new(fasta_file);

    let mut header: String = String::new();
    let mut aa_sequence = String::new();
    let trypsin: Trypsin = Trypsin::new(2, 6, 50);

    let mut peptides_from_fasta: Collection<Peptide> = Collection::new();

    for line in fasta_file.lines() {
        // trim
        let string_line = line.unwrap().as_mut_str().trim().to_owned()  ;
        if !string_line.starts_with(">") {
            aa_sequence.push_str(&string_line);
        } else {
            if header.len() > 0 {
                let mut protein: Protein = Protein::new(header.clone(), aa_sequence);
                // throw error if aa sequence has whitespaces
                trypsin.digest(&mut protein, &mut peptides_from_fasta);
                aa_sequence = String::new();
            }
            header = string_line;
        }
    }
    // process last protein
    let mut protein: Protein = Protein::new(header, aa_sequence);
    // throw error if aa sequence has whitespaces
    trypsin.digest(&mut protein, &mut peptides_from_fasta);


    // read peptides from tsv files which is digested by a reference program
    // change path so it fits tsv file
    path.set_file_name("digest_results_ju");
    path.set_extension("tsv");

    let tsv_file = File::open(path.to_str().unwrap()).expect("tsv file not found");
    let tsv_file = BufReader::new(tsv_file);

    let mut first_line: bool = true;    // first line contains headers
    let mut peptides_from_tsv: Collection<Peptide> = Collection::new();

    for line in tsv_file.lines() {
        // ignore first line
        if !first_line {
            let string_line = line.unwrap().as_mut_str().trim().to_owned();
            let splited_line: Vec<&str> = string_line.split('\t').collect();
            let peptide: Peptide = Peptide::new(String::from(splited_line[0]), String::from("Compare tsv file"), 0, 0);
            peptides_from_tsv.add(peptide);
        } else {
            first_line = false;
        }
    }

    // check if the length of the collection is the same
    println!("peptides_from_fasta_file => {}", peptides_from_fasta.len());
    println!("peptides_from_tsv_file => {}", peptides_from_tsv.len());
    assert_eq!(peptides_from_fasta.len(), peptides_from_tsv.len());
}
