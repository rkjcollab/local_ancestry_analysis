#!/bin/bash

# TO NOTE: update the paths below to set the path to the ADMIXTURE input
# directory (admix_input_dir), the RFMix results directory (rfmix_results_dir),
# and the output directory(out_dir).

# TODO: for now, must use SNPsPerWindow & Viterbi files from -co run. Could
# test making pseudo SNPsPerWindow file from -o run in future.

# Scripts assume RFMix output file names are as output by RFMix:
    # chr#_local_ancestry.0.Viterbi.txt
    # chr#_local_ancestry_batch1.0.SNPsPerWindow.txt

admix_dir="data/admixture_unimp"
rfmix_results_dir="data/rfmix_unimp/output_co"
out_dir="data/gwide_anc_unimp"

# Calculate local AFR ancestry
plink_file_name="${admix_dir}/input/initial_admixed"
nr_samples=`wc -l ${plink_file_name}.fam | xargs | cut -f1 -d' '`

Rscript code/calc_local_afr_ancestry/calc_local_afr_ancestry.R \
    $nr_samples $rfmix_results_dir $out_dir
