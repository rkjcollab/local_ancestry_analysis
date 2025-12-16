#!/bin/bash

chr=$1  # either pass single chromosome number or column for bed file with chromosomes
in_file=$2
out_dir=$3  # SDS added

# Get current code dir
code_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Puts into bed format for CrossMap
    # Need chrom, start, end
    # Getting error due to first row having "position", SDS adding code
    # to remove that

# input passed for single chromosome
if [[ $chr =~ ^[0-9]+$ ]]; then 
    chr="chr${chr}"
    cut -f1 $in_file > "$out_dir"/tmp_${chr}_c2.txt
    paste -d' ' "$out_dir"/tmp_${chr}_c2.txt "$out_dir"/tmp_${chr}_c2.txt "$out_dir"/tmp_${chr}_c2.txt \
        > "$out_dir"/tmp_c2_3_4.txt
    sed "s/^/$chr /" "$out_dir"/tmp_c2_3_4.txt | \
        tail -n +2 > "$out_dir"/tmp_${chr}_in.bed  # SDS added code here

    CrossMap bed ${code_dir}/hg19ToHg38.over.chain \
                "$out_dir"/tmp_${chr}_in.bed  \
                "$out_dir"/tmp_${chr}_out.bed

# input passed with chromosome column for bed file
else
    cut -f1 $in_file > "$out_dir"/tmp_c2.txt
    paste -d' ' "$chr" "$out_dir"/tmp_c2.txt "$out_dir"/tmp_c2.txt "$out_dir"/tmp_c2.txt \
        > "$out_dir"/tmp_c1_c2_3_4.txt

    # SDS added code here
    # Only use tail if grep finds alpha characters in first row
    if head -n 1 "${out_dir}"/tmp_c1_c2_3_4.txt | grep -q '[A-Za-z]'; then
        tail -n +2 "${out_dir}"/tmp_c1_c2_3_4.txt > "$out_dir"/tmp_in.bed
    else
        cat "${out_dir}"/tmp_c1_c2_3_4.txt > "$out_dir"/tmp_in.bed
    fi

    CrossMap bed ${code_dir}/hg19ToHg38.over.chain \
                "$out_dir"/tmp_in.bed  \
                "$out_dir"/tmp_out.bed

fi
