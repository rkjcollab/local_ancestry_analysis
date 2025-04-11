#!/bin/bash

# TO NOTE: First, set the paths to the data directory and the code directory, which
# will be bound to the container. Then update the input_file_name path to the same
# input used in the first step of the local_anc_afr_eur subfolder:
    # Unphased data: set input_file_name to plink_input in job_prep_input_unphased.batch,
        # does not matter which chromosome or if all chromosomes concatenated
    # Phased data: set input_file_name to input_vcf in job_prep_input_phased.batch,
        # does not matter which chromosome or if all chromosomes concatenated
# Finally, set the path to RFMix results split by chromosome (rfmix_results dir).
# TO NOTE: this script needs to be called from laaa_fine_mapping, not batch_files.

# Test WGS BARD
data_dir="/scratch/alpine/sslack@xsede.org/ortega/pipeline_test_data/bard_from_abhishek/rfmix"
code_dir="/projects/sslack@xsede.org/repos/local_ancestry_analysis/laaa_fine_mapping"
rfmix_results_dir="${data_dir}/imp_output_o"

# Test full imputed SARP/CSGA
# data_dir="/Users/slacksa/Library/CloudStorage/OneDrive-TheUniversityofColoradoDenver/Collabs/ortega/data/rfmix_wgs_new_2"
# code_dir="/Users/slacksa/repos/local_ancestry_analysis/laaa_fine_mapping"
# rfmix_results_dir="${data_dir}/output_o"

# Run script
cont_dir=$(dirname "${code_dir}")
apptainer exec --bind ${data_dir}:${data_dir} --bind ${code_dir}:${code_dir} \
    ${cont_dir}/local_ancestry_analysis.sif \
    Rscript ${code_dir}/calc_gwide_ancestry/calc_rfmix_gwide_ancestry.R \
        $rfmix_results_dir
