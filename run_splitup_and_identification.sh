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


shopt -s nullglob   # prevent expansion to *.mzML if not files are present, see: *https://www.cyberciti.biz/faq/bash-loop-over-file/
for spectrum_file in ${RESULTS_FOLDER}/*.mzML
do
    comet_result_dir=${RESULTS_FOLDER}/${spectrum_file}.d/
    echo "`date +%Y-%m-%d` start with ${spectrum_file}"
    /home/max-decoy/mac_decoy identification \
        --modification-file $MODIFICATION_FILE
        --spectrum-file $spectrum_file
        --max-number-of-variable-modification-per-peptide $MAX_NUMBER_OF_VARIABLE_MODIFICATION_PER_PEPTID
        --number-of-decoys $NUMBER_OF_DECOYS
        --lower-mass-tolerance $LOWER_MASS_TOLERANCE_IN_PPM
        --upper-mass-tolerance $UPPER_MASS_TOLERANCE_IN_PPM
        --thread-count $THREAD_COUNT
    /home/max-decoy/crux ./crux comet \
        --decoy_search 0 \
        --num_results 100 \
        --add_C_cysteine 0 \
        --num_output_lines 100 \
        --output-dir $comet_result_dir
        --output_txtfile 1 \
        --output_outfiles 1 \
        $spectrum_file \
        ${spectrum_file}.fasta
    echo "`date +%Y-%m-%d` with ${spectrum_file}"
    echo ""
    echo ""
done
