#!/bin/bash

# This script only calculates global ancestry from RFMix estimates, to be used as a 
# covariate in LAAA fine mapping. RFMix results for ALL chromosomes must be generated
# before this step is run. Can be run using RFMix results run with collapse flag -o
# or -co, but assumes RFMix was run in batches.

# TO NOTE: update the input_file_name path to the same input used in the first step
# of the local_anc_afr_eur subfolder:
    # Unphased data: set input_file_name to plink_input in job_prep_input_unphased.batch,
        # does not matter which chromosome or if all chromosomes concatenated
    # Phased data: set input_file_name to input_vcf in job_prep_input_phased.batch,
        # does not matter which chromosome or if all chromosomes concatenated
# Also set the path to RFMix results split by chromosome (rfmix_results dir).

# Scripts assume RFMix output file names are as output by RFMix:
    # chr#_local_ancestry.0.Viterbi.txt
    # chr#_local_ancestry.allelesRephased0.txt
    # chr#_local_ancestry_snps.txt
    # chr#_local_ancestry_samples.txt
    # chr#_local_ancestry_batch1.0.SNPsPerWindow.txt (if found in dir, means
        # RFMix was run with -co option)

# Test imputed
# input_file_name="${RKJCOLLAB}/Collabs/ortega/data/pipeline_test_data/imp/chr22_small.vcf"
# rfmix_results_dir="${RKJCOLLAB}/Collabs/ortega/data/pipeline_test_data/imp_output"

# Test chip
# input_file_name="${RKJCOLLAB}/Collabs/ortega/data/pipeline_test_data/chip/chr22_small"
# rfmix_results_dir="${RKJCOLLAB}/Collabs/ortega/data/pipeline_test_data/chip_output"

# Test WGS
input_file_name="${RKJCOLLAB}/Collabs/ortega/data/pipeline_test_data/local_anc_afr_eur/wgs/chr10_small.vcf"
rfmix_results_dir="${RKJCOLLAB}/Collabs/ortega/data/pipeline_test_data/local_anc_afr_eur/wgs_output"

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
