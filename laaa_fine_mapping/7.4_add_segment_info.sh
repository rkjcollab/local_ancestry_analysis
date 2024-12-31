#!/bin/bash

# TODO: as of 9/23/2024, not using this script

# TO NOTE: update the paths below to set ...

# Script gets admixture mapping regions in hg19 and converts them to hg38,
# to allow for annotation of LAAA results by admixture mapping segment.

# TODO: add grep command in liftover_to_hg38.sh
# TODO: confirm okay to assume order will not change
# TODO: remove region ext once can look up

# Set inputs
proj_dir="${RKJCOLLAB}/Collabs/ortega"
out_dir="${proj_dir}/data/pheno_analysis/output"
pheno=maxFEV1_FVC  # switch pheno to run for each region
chr=10  # switch chr to run for each region
# pheno=maxFVC
# chr=5
# chr=8
# chr=17

# Region file to be filled in with above phenotype and chromosome, made by admixture
# mapping summary report with one peak and its contiguous region
region="${proj_dir}/data/pheno_analysis/output/${pheno}_admixture_peak_contig_region_chr${chr}.txt"

# Get each region start and stop for liftover
cut -f3 $region > "${out_dir}/tmp_range.txt"
cut -d '-' -f1 "${out_dir}/tmp_range.txt" > "${out_dir}/tmp_range_start.txt"
cut -d '-' -f2 "${out_dir}/tmp_range.txt" > "${out_dir}/tmp_range_stop.txt"

# Liftover
# Output automatically formatted as tmp_chr${chr}_out.bed
bash make_dose_frames/liftover_to_hg38.sh \
    "$chr" \
    "${out_dir}/tmp_range_start.txt" \
    "$out_dir"
cut -f2 "${out_dir}/tmp_chr${chr}_out.bed" > "${out_dir}/tmp_range_start_hg38.txt"
echo "BEGIN.HG38" > "${out_dir}/tmp_range_start_head.txt"
cat "${out_dir}/tmp_range_start_head.txt" \
    "${out_dir}/tmp_range_start_hg38.txt" > \
    "${out_dir}/tmp_range_start_hg38_head.txt"

bash make_dose_frames/liftover_to_hg38.sh \
    "$chr" \
    "${out_dir}/tmp_range_stop.txt" \
    "$out_dir"
cut -f2 "${out_dir}/tmp_chr${chr}_out.bed" > "${out_dir}/tmp_range_stop_hg38.txt"
echo "END.HG38" > "${out_dir}/tmp_range_stop_head.txt"
cat "${out_dir}/tmp_range_stop_head.txt" \
    "${out_dir}/tmp_range_stop_hg38.txt" > \
    "${out_dir}/tmp_range_stop_hg38_head.txt"

# Combine and write out
# TODO: edit region format here!
region_out="${proj_dir}/data/pheno_analysis/output/${pheno}_admixture_peak_contig_region_chr${chr}_hg38.txt"
paste -d'\t' "$region" "${out_dir}/tmp_range_start_hg38_head.txt" \
    "${out_dir}/tmp_range_stop_hg38_head.txt" > \
    "$region_out"

# Clean up
rm ${out_dir}/tmp_*
