#!/bin/bash

# TO NOTE: update the paths below to set ...

# Regions are contiguous segments with P < 1E-3 around significant peaks.

# Set inputs
proj_dir="${RKJCOLLAB}/Collabs/ortega"
out_dir_prefix="${proj_dir}/data/laaa_wgs_new_2/output"
pheno=maxFEV1_FVC  # switch pheno to run for each region
chr=10  # switch chr to run for each region
# pheno=maxFVC
# chr=5
# chr=8
# chr=17

# Dose frame to be filled in with above phenotype and chromosome
dose_dir="${proj_dir}/data/laaa_wgs_new_2/dose_frames/region_hg38_${pheno}_chr${chr}"

# Pheno file to be filled in with above phenotype
pheno_file="${proj_dir}/data/pheno/SARP123_CSGA_348_laaa_${pheno}_pheno.txt"
pheno_id_col_name="TopMed_ID"  # should match format in RFMix output

# Run models for region
mkdir "$out_dir_prefix"
out_dir_base=$(basename "$dose_dir")
out_dir="${out_dir_prefix}/${out_dir_base}"
mkdir $out_dir

Rscript run_models/run_models.R \
    $pheno $chr $dose_dir $out_dir $pheno_file $pheno_id_col_name
