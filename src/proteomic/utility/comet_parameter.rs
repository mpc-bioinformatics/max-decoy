use std::collections::HashMap;
use std::path::Path;

use proteomic::models::amino_acids::modification::Modification;

const COMET_PARAMS_BEGIN: &'static str = "
# comet_version 2018.01 rev. 4
# Comet MS/MS search engine parameters file.
# Everything following the '#' symbol is treated as a comment.

decoy_search = 0
peff_format = 0
peff_obo =

num_threads = 0

peptide_mass_units = 2
mass_type_parent = 1
mass_type_fragment = 1
precursor_tolerance_type = 1
isotope_error = 3

search_enzyme_number = 1
num_enzyme_termini = 2
allowed_missed_cleavage = 2

max_variable_mods_in_peptide = 5
require_variable_mod = 0

theoretical_fragment_ions = 1
use_A_ions = 0
use_B_ions = 1
use_C_ions = 0
use_X_ions = 0
use_Y_ions = 1
use_Z_ions = 0
use_NL_ions = 0

output_sqtstream = 0
output_sqtfile = 0
output_txtfile = 1
output_pepxmlfile = 0
output_percolatorfile = 0
print_expect_score = 1
show_fragment_ions = 0

sample_enzyme_number = 1

scan_range = 0 0
precursor_charge = 0 0
override_charge = 0
ms_level = 2
activation_method = ALL

digest_mass_range = 600.0 5000.0
skip_researching = 1
max_fragment_charge = 3
max_precursor_charge = 6
nucleotide_reading_frame = 0
clip_nterm_methionine = 0
spectrum_batch_size = 0
decoy_prefix = DECOY_
equal_I_and_L = 1
output_suffix =
mass_offsets =

minimum_peaks = 10
minimum_intensity = 0
remove_precursor_peak = 0
remove_precursor_tolerance = 1.5
clear_mz_range = 0.0 0.0

add_Cterm_peptide = 0.0
add_Nterm_peptide = 0.0
add_Cterm_protein = 0.0
add_Nterm_protein = 0.0

fragment_bin_offset = 0

";

const COMET_PARAMS_END: &'static str = "
[COMET_ENZYME_INFO]
0.  No_enzyme              0      -           -
1.  Trypsin                1      KR          P
2.  Trypsin/P              1      KR          -
3.  Lys_C                  1      K           P
4.  Lys_N                  0      K           -
5.  Arg_C                  1      R           P
6.  Asp_N                  0      D           -
7.  CNBr                   1      M           -
8.  Glu_C                  1      DE          P
9.  PepsinA                1      FL          P
10. Chymotrypsin           1      FWYL        P
";

pub fn new(fix_modifications_map: &HashMap<char, Modification>, variable_modifications_map: &HashMap<char, Modification>, fasta_file_path: &Path, number_of_target_and_decoys: usize, max_number_of_variable_modification_per_peptide: u8, fragmentation_tolerance: f64, lower_precursor_tolerance: i64, upper_precursor_tolerance: i64) -> String {
    let mut params: String = COMET_PARAMS_BEGIN.to_owned();
    // Comet has only one parameter for upper and lower precursor tolerance which is peptide_mass_tolerance. So we use the greatest of them.
    params.push_str(format!("peptide_mass_tolerance = {:.4}\n", std::cmp::max(lower_precursor_tolerance, upper_precursor_tolerance) as f64).as_str());
    params.push_str(format!("fragment_bin_tol = {}\n", fragmentation_tolerance).as_str());
    params.push_str(format!("num_results = {}\n", number_of_target_and_decoys).as_str());
    params.push_str(format!("num_output_lines = {}\n", number_of_target_and_decoys).as_str());
    let fasta_file_path_as_str = match fasta_file_path.as_os_str().to_str() {
        Some(path) => path,
        None => panic!("proteomic::utility::comet_parameter::new(): No FASTA-file path")
    };
    params.push_str(format!("database_name = {}\n", fasta_file_path_as_str).as_str());
    for (_, modification) in fix_modifications_map {
        params.push_str(modification.to_comet_static_modification_param().as_str());
        params.push_str("\n");
    }
    if !fix_modifications_map.contains_key(&'J') { params.push_str("add_J_user_amino_acid = 113.08406\n"); }
    let mut modification_number: u8 = 1;
    for (_, modification) in variable_modifications_map {
        params.push_str(modification.to_comet_variable_modification_param(modification_number, max_number_of_variable_modification_per_peptide).as_str());
        params.push_str("\n");
        if modification_number == 9 { break; }
        modification_number += 1;
    }
    params.push_str(COMET_PARAMS_END);
    return params;
}


