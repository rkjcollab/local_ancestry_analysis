#!/bin/bash

# Script preps input for local ancestry pipeline, currently
# written to start from a single chromosome input in hg38.

# Get args
vcf_input=$1
out_dir=$2
code_dir=$3

# Set PLINK output
plink_output_temp=$(basename "$vcf_input")
plink_output_temp2="${plink_output_temp%.*}"
plink_output="${plink_output_temp2%.*}"  # repeat for .vcf.gz

# Prep input
# TODO: need to add PLINK2 to container, currently this is run in
# batch script
# plink2 \
# 	--vcf "$vcf_input" \
# 	--make-bed  --keep-allele-order \
# 	--set-all-var-ids @:#:\$r:\$a \
# 	--new-id-max-allele-len 100 \
# 	--out "${out_dir}/${plink_output}"

# Prep input
# Remove duplicate positions, keeping first (matches topmed_imputation repo)
plink \
	--bfile "${out_dir}/${plink_output}" \
	--list-duplicate-vars 'ids-only' 'suppress-first'\
	--out "${out_dir}/tmp_${plink_output}_dup_list"

plink \
	--bfile "${out_dir}/${plink_output}" \
	--make-bed --keep-allele-order \
	--exclude "${out_dir}/tmp_${plink_output}_dup_list.dupvar" \
	--out "${out_dir}/tmp_${plink_output}_dedup"

# Liftover, based on code from topmed_imputation repo
# Create bed file to crossover from hg38 to hg19
cat "${out_dir}/tmp_${plink_output}_dedup.bim" | cut -f1 | sed 's/^/chr/' > "${out_dir}/tmp_c1.txt"
cat "${out_dir}/tmp_${plink_output}_dedup.bim" | cut -f4 > "${out_dir}/tmp_c2.txt"
cat "${out_dir}/tmp_${plink_output}_dedup.bim" | cut -f4 > "${out_dir}/tmp_c3.txt"
cat "${out_dir}/tmp_${plink_output}_dedup.bim" | cut -f2 > "${out_dir}/tmp_c4.txt"
paste  "${out_dir}/tmp_c1.txt" \
       "${out_dir}/tmp_c2.txt" \
       "${out_dir}/tmp_c3.txt" \
       "${out_dir}/tmp_c4.txt" \
       > "${out_dir}/tmp_in.bed"

# Do crossover
CrossMap.py bed ${code_dir}/shapeit_formatting_scripts/hg38ToHg19.over.chain \
	"${out_dir}/tmp_in.bed"  \
	"${out_dir}/tmp_out.bed"

# Extract only those SNPs that were successfully cross-overed
cut -f4 "${out_dir}/tmp_out.bed" > "${out_dir}/tmp_snp_keep.txt"
plink \
	--bfile "${out_dir}/tmp_${plink_output}_dedup" \
	--extract "${out_dir}/tmp_snp_keep.txt" \
	--keep-allele-order \
	--make-bed --out "${out_dir}/${plink_output}_dedup_hg19"

# Update bim file positions
Rscript --vanilla ${code_dir}/update_pos.R \
	"${out_dir}/tmp_out.bed" "${out_dir}/${plink_output}_dedup_hg19.bim"

# Clean up
rm ${out_dir}/tmp_*
