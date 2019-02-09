#!/bin/sh

SPECTRA_FILE=/home/max-decoy/spectra.mzML
MODIFICATION_FILE=/home/max-decoy/modifications.csv
RESULTS_FOLDER=/home/max-decoy/results/

curl -o $SPECTRA_FILE $SPECTRA_FILE_URL
curl -o $MODIFICATION_FILE $MODIFICATION_FILE_URL

/home/max-decoy/mac_decoy spectrum-splitup \
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
    echo "`date +%Y-%m-%d` start with ${spectrum_file_path}"
    /home/max-decoy/mac_decoy identification \
        --modification-file $MODIFICATION_FILE \
        --spectrum-file $spectrum_file_path \
        --max-number-of-variable-modification-per-peptide $MAX_NUMBER_OF_VARIABLE_MODIFICATION_PER_PEPTID \
        --number-of-decoys $NUMBER_OF_DECOYS \
        --lower-mass-tolerance $LOWER_MASS_TOLERANCE_IN_PPM \
        --upper-mass-tolerance $UPPER_MASS_TOLERANCE_IN_PPM \
        --thread-count $THREAD_COUNT
    /home/max-decoy/crux comet \
        --decoy_search 0 \
        --num_results 100 \
        --add_C_cysteine 0 \
        --num_output_lines 100 \
        --output-dir $comet_result_dir \
        --output_txtfile 1 \
        --output_outfiles 1 \
        $spectrum_file_path \
        ${RESULTS_FOLDER}/${spectrum_filename}.fasta
    echo "`date +%Y-%m-%d` with ${spectrum_file_path}"
    echo ""
    echo ""
done
