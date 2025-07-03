#!/bin/bash

# Runs LAAA models for each region present in dose frame directory. Only
# runs models for hg38 dose frames.

out_dir_prefix="$1"
dose_dir="$2"
pheno_file="$3"
pheno_id_col_name="$4" # should match format in RFMix output
cov_list="$5"  # should be passed as "cov1,cov2" column names in pheno file

# Loop over dose frames and run LAAA for each
dose_dir_list=$(ls $dose_dir | grep "region_hg38")
for dose in $dose_dir_list; do
    echo "Running LAAA for ${dose}."

    [[ $dose =~ region_hg38_(.*)_chr(.*) ]];
    pheno=${BASH_REMATCH[1]}
    chr=${BASH_REMATCH[2]}

    # Pheno file to be filled in with above phenotype
    pheno_file_pheno=$(eval "echo $pheno_file")

    # Run models for region
    mkdir "$out_dir_prefix"
    out_dir="${out_dir_prefix}/${dose}"
    mkdir $out_dir

    echo $chr

    Rscript run_models/run_models.R \
        $pheno $chr ${dose_dir}/${dose} $out_dir $pheno_file_pheno $pheno_id_col_name $cov_list

done
