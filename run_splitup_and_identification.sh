#!/bin/bash

# This script uses the tide-version of comet for CLI-support.
# Adjust the tide-comet-params to your needs

# executables
MAX_DECOY=./max_decoy
TIDE=./tide

# folder settings
SPECTRA_FILE=./spectra.mzML
MODIFICATION_FILE=./test_modifications.csv
RESULTS_FOLDER=./results
SPECTRUM_SUFFIX=""

# identification settings
MAX_NUMBER_OF_VARIABLE_MODIFICATION_PER_PEPTID=2
NUMBER_OF_DECOYS=5000
LOWER_MASS_TOLERANCE_IN_PPM=5
UPPER_MASS_TOLERANCE_IN_PPM=5
THREAD_COUNT=30
TIMEOUT=60

$MAX_DECOY spectrum-splitup \
    --mz-ml-file $SPECTRA_FILE \
    --destination-folder $RESULTS_FOLDER \
    --file-suffix $SPECTRUM_SUFFIX


shopt -s nullglob   # prevent expansion to *.mzML if no files are present, see: *https://www.cyberciti.biz/faq/bash-loop-over-file/
for spectrum_file_path in ${RESULTS_FOLDER}/*.mzML
do
    spectrum_filename=$(basename -- "$spectrum_file_path")
    spectrum_filename="${spectrum_filename%.*}"
    comet_result_dir=${RESULTS_FOLDER}/${spectrum_filename}_results/
    mkdir $comet_result_dir
    echo "`date +%Y-%m-%d\ %H:%M:%S` start with ${spectrum_file_path}"
    $MAX_DECOY identification \
        --modification-file $MODIFICATION_FILE \
        --spectrum-file $spectrum_file_path \
        --max-number-of-variable-modification-per-peptide $MAX_NUMBER_OF_VARIABLE_MODIFICATION_PER_PEPTID \
        --number-of-decoys $NUMBER_OF_DECOYS \
        --lower-mass-tolerance $LOWER_MASS_TOLERANCE_IN_PPM \
        --upper-mass-tolerance $UPPER_MASS_TOLERANCE_IN_PPM \
        --thread-count $THREAD_COUNT \
    	--max-time-for-decoy-generation $TIMEOUT
    $TIDE comet \
	--fragment_bin_tol 0.02 \
	--fragment_bin_offset 0 \
        --decoy_search 0 \
        --num_results 100 \
        --add_C_cysteine 57.021464 \
        --add_J_user_amino_acid 113.08406 \
	    --variable_mod01 = 15.994915 M 0 3 -1 0 0 \
        --num_output_lines 100 \
        --output-dir $comet_result_dir \
        --output_txtfile 1 \
        --output_outfiles 1 \
	    --overwrite T \
        $spectrum_file_path \
        ${RESULTS_FOLDER}/${spectrum_filename}.fasta
    cp ${comet_result_dir}comet.target.txt ${RESULTS_FOLDER}/${spectrum_filename}.result.tsv
    rm -r $comet_result_dir
    echo "`date +%Y-%m-%d\ %H:%M:%S` with ${spectrum_file_path}"
    echo ""
    echo ""
done

