#!/bin/bash

chr=$1  # either pass single chromosome number or column for bed file with chromosomes
in_file=$2
out_dir=$3  # SDS added

# Puts into bed format for CrossMap
    # Need chrom, start, end
    # Getting error due to first row having "position", SDS adding code
    # to remove that

# TODO: Need to update to only use tail if grep finds alpha characters in first row

# input passed for single chromosome
if [[ $chr =~ ^[0-9]+$ ]]; then 
    chr="chr${chr}"
    cut -f1 $in_file > "$out_dir"/tmp_${chr}_c2.txt
    paste -d' ' "$out_dir"/tmp_${chr}_c2.txt "$out_dir"/tmp_${chr}_c2.txt "$out_dir"/tmp_${chr}_c2.txt \
        > "$out_dir"/tmp_c2_3_4.txt
    sed "s/^/$chr /" "$out_dir"/tmp_c2_3_4.txt | \
        tail -n +2 > "$out_dir"/tmp_${chr}_in.bed  # SDS added code here

    CrossMap bed make_dose_frames/hg19ToHg38.over.chain \
                "$out_dir"/tmp_${chr}_in.bed  \
                "$out_dir"/tmp_${chr}_out.bed

# input passed with chromosome column for bed file
else
    cut -f1 $in_file > "$out_dir"/tmp_c2.txt
    paste -d' ' "$chr" "$out_dir"/tmp_c2.txt "$out_dir"/tmp_c2.txt "$out_dir"/tmp_c2.txt \
        > "$out_dir"/tmp_c1_c2_3_4.txt
    tail -n +2 "${out_dir}"/tmp_c1_c2_3_4.txt > "$out_dir"/tmp_in.bed  # SDS added code here

    CrossMap bed make_dose_frames/hg19ToHg38.over.chain \
                "$out_dir"/tmp_in.bed  \
                "$out_dir"/tmp_out.bed

fi
