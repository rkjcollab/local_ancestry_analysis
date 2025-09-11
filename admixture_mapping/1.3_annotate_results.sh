#!/bin/bash

# TODO: add description here matching rest of pipeline.

# TO NOTE: activate crossmap-bcftools-osx64 before running.

# Set inputs
plink_file_name=$1
rfmix_dir=$2
result_dir=$3
pheno_dir=$4

# Get list of phenotype files input into model, use extracted
# phenotype to get results file name
pheno_file_list=$(ls $pheno_dir | grep -E "pheno(_male|_female)?\.txt")
for pheno_file in $pheno_file_list; do
    [[ $pheno_file =~ SARP123_CSGA_[0-9]+_admix_map_(.*)_pheno(_male|_female)?\.txt ]];
    pheno=${BASH_REMATCH[1]}
    sex=${BASH_REMATCH[2]}  # may be empty

    # Get assoc results for that pheno
    assoc_file="SARP123_CSGA_${pheno}_admix_map${sex}.assoc.linear"

    Rscript admix_mapping/annotate_results.R \
        "${plink_file_name}.fam" \
        "$rfmix_dir" \
        "${pheno_dir}/${pheno_file}" \
        "${result_dir}/${assoc_file}"

    # Make version of each full results file in hg38
    # Get each region start and stop for liftover
    cut -f3 "${result_dir}/${assoc_file}_annot.txt" > "${result_dir}/tmp_range.txt"
    cut -d '-' -f1 "${result_dir}/tmp_range.txt" > "${result_dir}/tmp_range_start.txt"
    cut -d '-' -f2 "${result_dir}/tmp_range.txt" > "${result_dir}/tmp_range_stop.txt"
    cut -f1 "${result_dir}/${assoc_file}_annot.txt" > "${result_dir}/tmp_chr.txt"

    # Liftover
    # Output automatically formatted as tmp_chr${chr}_out.bed

    chr="${result_dir}/tmp_chr.txt"  # pass bed file with chromosome positions
    bash make_dose_frames/liftover_to_hg38.sh \
        "$chr" \
        "${result_dir}/tmp_range_start.txt" \
        "$result_dir"
    cut -f2 "${result_dir}/tmp_out.bed" > "${result_dir}/tmp_range_start_hg38.txt"
    echo "BEGIN.HG38" > "${result_dir}/tmp_range_start_head.txt"
    cat "${result_dir}/tmp_range_start_head.txt" \
        "${result_dir}/tmp_range_start_hg38.txt" > \
        "${result_dir}/tmp_range_start_hg38_head.txt"

    bash make_dose_frames/liftover_to_hg38.sh \
        "$chr" \
        "${result_dir}/tmp_range_stop.txt" \
        "$result_dir"
    cut -f2 "${result_dir}/tmp_out.bed" > "${result_dir}/tmp_range_stop_hg38.txt"
    echo "END.HG38" > "${result_dir}/tmp_range_stop_head.txt"
    cat "${result_dir}/tmp_range_stop_head.txt" \
        "${result_dir}/tmp_range_stop_hg38.txt" > \
        "${result_dir}/tmp_range_stop_hg38_head.txt"

    # Combine and write out
    # TODO: edit region format here!
    assoc_file_out="${result_dir}/${assoc_file}_annot_hg38.txt"
    paste -d'\t' "${result_dir}/${assoc_file}_annot.txt" "${result_dir}/tmp_range_start_hg38_head.txt" \
        "${result_dir}/tmp_range_stop_hg38_head.txt" > \
        "$assoc_file_out"

    # Clean up
    rm ${result_dir}/tmp_*

done
