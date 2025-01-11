#!/bin/bash

# TO NOTE: First, set the paths to the data directory and the code directory, which
# will be bound to the container. Then, set inputs for making dose frames:
    # the output directory prefix (out_dir_prefix), the path to RFMix results split
    # by chromosome (rfmix_results dir),and the path to the directory containing
    # admixture mapping peak regions (admix_dir) in hg19 (written out by admixture
    # mapping summary report).
# Finally, set inputs for running LAAA models:
    # the phenotype file for analysis (pheno_file), the name of the column with
    # IDs in the phenotype file in a format that matches RFMix output
    # (pheno_id_col_name), and the names of the columns in the phenotype file
    # to use as covariates (cov_list).

# Script outputs dose frames in hg19 and hg38, and LAAA is run in hg38.

# Set inputs for making dose frames
data_dir="/Users/slacksa/Library/CloudStorage/OneDrive-TheUniversityofColoradoDenver/Collabs/ortega/data/pipeline_test_data"
code_dir="/Users/slacksa/repos/local_ancestry_analysis/laaa_fine_mapping"
out_dir_prefix="${data_dir}/laaa_fine_mapping/laaa_wgs"
    # Script makes sub-directories "dose_frames" and "output" under out_dir_prefix
rfmix_results_dir="${data_dir}/local_anc_afr_eur/wgs_output"
admix_dir="${data_dir}/laaa_fine_mapping/pheno_analysis/output"

# Set inputs for running LAAA model
pheno_file='laaa_fine_mapping/laaa_wgs/pheno/SARP123_CSGA_348_laaa_${pheno}_pheno.txt'
    # single quotes are required here for substitution of and $pheno below
    # pheno string in name must match column name with phenotype values
pheno_id_col_name="TopMed_ID"  # should match format in RFMix output
cov_list="group,age,sex,ht,bmi,RFMIX_GW_AFR"  # should be passed as "cov1,cov2" column names in pheno file

# Run step 2.2 script to make dose frames
cont_dir=$(dirname "${code_dir}")
apptainer exec --bind ${data_dir}:${data_dir} --bind ${code_dir}:${code_dir} \
    ${cont_dir}/local_ancestry_analysis.sif \
    bash ${code_dir}/2.2_create_dose_frames.sh \
    "${out_dir_prefix}/dose_frames" \
    "$rfmix_results_dir" \
    "$admix_dir"

# Run step 2.3 script to run LAAA model
out_dir="${out_dir_prefix}/output"
dose_dir="${out_dir_prefix}/dose_frames"
apptainer exec --bind ${data_dir}:${data_dir} --bind ${code_dir}:${code_dir} \
    ${cont_dir}/local_ancestry_analysis.sif \
    bash ${code_dir}/2.3_run_models.sh \
    "${out_dir_prefix}/output" \
    "${out_dir_prefix}/dose_frames" \
    "${data_dir}/${pheno_file}" \
    "$pheno_id_col_name" \
    "$cov_list"
