#!/bin/bash

chr=$1

#Get lifted over coordinates
cat get_snps_to_include.R | R --vanilla --args $chr
CrossMap.py bed hg38ToHg19.over.chain \
            tmp_chr${chr}_in.bed  \
            tmp_chr${chr}_out.bed
cat rm_dupl_pos.R | R --vanilla --args $chr

#Process files
python create_shapeit_formatted_output.py $chr

#Cleanup
rm tmp_chr${chr}_in.bed
rm tmp_chr${chr}_out.bed
