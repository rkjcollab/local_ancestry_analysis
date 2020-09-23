#!/bin/bash

dose_file=$1
info_file=$2
min_rsq=$3
min_maf=$4

#Get chromsome nr
chr=`zcat $dose_file | tail -1 |  awk '{print $1}' | sed 's/chr//'`

#Get lifted over coordinates
cat /home/shapeit_formatting_scripts/get_imputed_snps_to_include.R | R --vanilla --args $info_file tmp_chr${chr}_in.bed $min_rsq $min_maf 
CrossMap.py bed /home/shapeit_formatting_scripts/hg38ToHg19.over.chain \
            tmp_chr${chr}_in.bed  \
            tmp_chr${chr}_out.bed
cat /home/shapeit_formatting_scripts/rm_crossmap_dupl_pos_order_by_hg19.R | R --vanilla --args $chr tmp_chr${chr}_out.bed

#Process files
python /home/shapeit_formatting_scripts/create_shapeit_formatted_output.py $chr $dose_file tmp_chr${chr}_out.bed
