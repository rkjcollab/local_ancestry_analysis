#!/bin/bash

# TO NOTE: update the paths below to set the project directory (proj_dir),
# the output directory (out_dir_prefix), the directory containing
# the dose frame region folders (dose_dir), and the phenotype file for
# analysis (pheno_file). 

# TO NOTE: update the paths below to set the chromosome (chr), windows
# (windows, list of "start stop" in hg19), phenotype for analysis (pheno)
# the input data directory (data_dir), and the output directory prefix
# (out_dir_prefix). Outputs dose frames in hg19 and hg38. Also set path
# admix_dir to directory containing admixture mapping peak regions in
# hg19 (written out by summary report).

# Regions are contiguous segments with P < 1E-3 around significant peaks.

# Set inputs
# Test WGS
proj_dir="${RKJCOLLAB}/Collabs/ortega"
out_dir_prefix="${proj_dir}/data/pipeline_test_data/laaa_fine_mapping/laaa_wgs/output"
dose_dir="${proj_dir}/data/pipeline_test_data/laaa_fine_mapping/laaa_wgs/dose_frames"
pheno_file='${proj_dir}/data/pheno/SARP123_CSGA_348_laaa_${pheno}_pheno.txt'
    # single quotes are required here for substitution of $proj_dir and $pheno below
    # pheno string in name must match column name with phenotype values
pheno_id_col_name="TopMed_ID"  # should match format in RFMix output
cov_list="group,age,sex,ht,bmi,RFMIX_GW_AFR"  # should be passed as "cov1,cov2" column names in pheno file

# Loop over dose frames and run LAAA for each
dose_dir_list=$(ls $dose_dir | grep "region_hg19")
for dose in $dose_dir_list; do
    echo "Running LAAA for ${dose}."

    [[ $dose =~ region_hg19_(.*)_chr(.*) ]];
    pheno=${BASH_REMATCH[1]}
    chr=${BASH_REMATCH[2]}

    # Pheno file to be filled in with above phenotype
    pheno_file_pheno=$(eval "echo $pheno_file")

    # Run models for region
    mkdir "$out_dir_prefix"
    out_dir_base=$(basename "$dose_dir")
    out_dir="${out_dir_prefix}/${out_dir_base}"
    mkdir $out_dir

    Rscript run_models/run_models.R \
        $pheno $chr ${dose_dir}/${dose} $out_dir $pheno_file_pheno $pheno_id_col_name $cov_list

done

