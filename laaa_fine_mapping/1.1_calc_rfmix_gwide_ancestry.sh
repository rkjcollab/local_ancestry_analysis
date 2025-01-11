#!/bin/bash

# This script only calculates global ancestry from RFMix estimates, to be used as a 
# covariate in LAAA fine mapping. RFMix results for ALL chromosomes must be generated
# before this step is run. Can be run using RFMix results run with collapse flag -o
# or -co. Assumes RFMix was run in batches.

# Scripts assume RFMix output file names are as output by RFMix:
    # chr#_local_ancestry.0.Viterbi.txt
    # chr#_local_ancestry.allelesRephased0.txt
    # chr#_local_ancestry_snps.txt
    # chr#_local_ancestry_samples.txt
    # chr#_local_ancestry_batch1.0.SNPsPerWindow.txt (if found in dir, means
        # RFMix was run with -co option)

input_file_name="$1"
rfmix_results_dir="$2"

# Get list of sample IDs from VCF or PLINK file
if [[ "$input_file_name" == *.vcf* ]]; then
    input_file_prefix="${input_file_name%.gz}"
    input_file_prefix="${input_file_prefix%.vcf}"
    bcftools query -l "$input_file_name" > \
        "${input_file_prefix}_sample_list.txt"
else
    awk '{print $2}' "${input_file_name}.fam" > "${input_file_name}_sample_list.txt"
fi

# Calculate RFMIX genome-wide ancestry
Rscript calc_gwide_ancestry/calc_rfmix_gwide_ancestry.R \
    "${input_file_prefix}_sample_list.txt" $rfmix_results_dir
