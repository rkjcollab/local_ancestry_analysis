#!/bin/bash

# TO NOTE: First, set the paths to the data directory and the code directory, which
# will be bound to the container. Then update the input_file_name path to the same
# input used in the first step of the local_anc_afr_eur subfolder:
    # Unphased data: set input_file_name to plink_input in job_prep_input_unphased.batch,
        # does not matter which chromosome or if all chromosomes concatenated
    # Phased data: set input_file_name to input_vcf in job_prep_input_phased.batch,
        # does not matter which chromosome or if all chromosomes concatenated
# Finally, set the path to RFMix results split by chromosome (rfmix_results dir).

# Test imputed
# input_file_name="${RKJCOLLAB}/Collabs/ortega/data/pipeline_test_data/imp/chr22_small.vcf"
# rfmix_results_dir="${RKJCOLLAB}/Collabs/ortega/data/pipeline_test_data/imp_output"

# Test chip
# input_file_name="${RKJCOLLAB}/Collabs/ortega/data/pipeline_test_data/chip/chr22_small"
# rfmix_results_dir="${RKJCOLLAB}/Collabs/ortega/data/pipeline_test_data/chip_output"

# Test WGS
# input_file_name="${RKJCOLLAB}/Collabs/ortega/data/pipeline_test_data/local_anc_afr_eur/wgs/chr10_small.vcf"
# rfmix_results_dir="${RKJCOLLAB}/Collabs/ortega/data/pipeline_test_data/local_anc_afr_eur/wgs_output"

# Test BARD
data_dir="/scratch/alpine/sslack@xsede.org/ortega/pipeline_test_data/bard_from_abhishek"
code_dir="/projects/sslack@xsede.org/repos/local_ancestry_analysis/laaa_fine_mapping"
input_file_name="${data_dir}/imp/chr1.dose.first10krows.vcf"
rfmix_results_dir="${data_dir}/imp_output"

# Run script
cont_dir=$(dirname "${code_dir}")
apptainer exec --bind ${data_dir}:${data_dir} --bind ${code_dir}:${code_dir} \
		${cont_dir}/local_ancestry_analysis.sif \
 		bash ${code_dir}/1.1_calc_rfmix_gwide_ancestry.sh \
            "$input_file_name" "$rfmix_results_dir"