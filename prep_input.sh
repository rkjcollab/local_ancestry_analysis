#!/bin/bash

# Script preps input for local ancestry pipeline, currently
# written to start from a single chromosome input in hg38.

# Get args
vcf_input=$1
plink_output=$2
out_dir=$3
code_dir=$4

# Prep input
# TODO: this PLINK2 not present in container!
# plink2 \
# 	--vcf "$vcf_input" \
# 	--make-bed  --keep-allele-order \
# 	--set-all-var-ids @:#:\$r:\$a \
# 	--new-id-max-allele-len 67 \
# 	--out "$plink_output"

# TODO: for now, arbitrarily keeping first of duplicate positions
plink \
	--bfile "$plink_output" \
	--list-duplicate-vars 'ids-only' 'suppress-first'\
	--out "${plink_output}_dup_list"

plink \
	--bfile "$plink_output" \
	--make-bed --keep-allele-order \
	--exclude "${plink_output}_dup_list.dupvar" \
	--out "${plink_output}_dedup"

# Based on code from topmed_imputation repo
#Create bed file to crossover from hg19 to hg38 
cat ${plink_output}_dedup.bim | cut -f1 | sed 's/^/chr/' > ${out_dir}/tmp_c1.txt
cat ${plink_output}_dedup.bim | cut -f4 > ${out_dir}/tmp_c2.txt
cat ${plink_output}_dedup.bim | cut -f4 > ${out_dir}/tmp_c3.txt
cat ${plink_output}_dedup.bim | cut -f2 > ${out_dir}/tmp_c4.txt
paste  ${out_dir}/tmp_c1.txt \
       ${out_dir}/tmp_c2.txt \
       ${out_dir}/tmp_c3.txt \
       ${out_dir}/tmp_c4.txt \
       >  ${out_dir}/tmp_in.bed

#Do crossover
CrossMap.py bed ${code_dir}/shapeit_formatting_scripts/hg38ToHg19.over.chain \
	${out_dir}/tmp_in.bed  \
	${out_dir}/tmp_out.bed

#Extract only those SNPs that were successfully cross-overed
cut -f4 ${out_dir}/tmp_out.bed > ${out_dir}/tmp_snp_keep.txt
plink \
	--bfile ${plink_output}_dedup \
	--extract ${out_dir}/tmp_snp_keep.txt \
	--keep-allele-order \
	--make-bed --out ${plink_output}_dedup_hg19

#Update bim file positions
Rscript --vanilla ${code_dir}/update_pos.R \
	${out_dir}/tmp_out.bed ${plink_output}_dedup_hg19.bim
