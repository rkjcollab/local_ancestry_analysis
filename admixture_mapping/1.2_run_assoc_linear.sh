#!/bin/bash

# Based on: Victor_Ortega_analyses_MPB/SARP_CSGA_merge/pheno_analysis/scripts/
#               run_assoc_linear_all_pheno_updated_COV.sh

# Set inputs
plink_input=$1
covar=$2
pheno_dir=$3
out_dir=$4

#TODO: need to remove hardcoded pheno name!

# Get list of phenotype files want to run analysis for
pheno_file_list=$(ls $pheno_dir | grep -E "pheno(_male|_female)?\.txt")
for pheno_file in $pheno_file_list; do
    [[ $pheno_file =~ SARP123_CSGA_[0-9]+_admix_map_(.*)_pheno(_male|_female)?\.txt ]];
    pheno=${BASH_REMATCH[1]}
    sex=${BASH_REMATCH[2]}  # may be empty

    out_file="${out_dir}/SARP123_CSGA_${pheno}_admix_map${sex}"

    plink --bfile "$plink_input" \
        --pheno "${pheno_dir}/${pheno_file}" \
        --covar "$covar" \
        --linear --out "$out_file"

done
