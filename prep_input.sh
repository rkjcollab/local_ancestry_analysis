#!/bin/bash

# Script preps input for local ancestry pipeline, currently
# written to start from a single chromosome PLINK input in either
# hg19 or hg38.

# Get args
plink_input=$1
out_dir=$2
code_dir=$3
build=$4

# Set PLINK output
plink_output=$(basename "$plink_input")

# Prep input
# Remove duplicate positions, keeping first
plink \
	--bfile "${plink_input}" \
	--keep-allele-order \
	--list-duplicate-vars 'ids-only' 'suppress-first'\
	--out "${out_dir}/tmp_${plink_output}_dup_list"

plink \
	--bfile "${plink_input}" \
	--make-bed --keep-allele-order \
	--exclude "${out_dir}/tmp_${plink_output}_dup_list.dupvar" \
	--out "${out_dir}/tmp_${plink_output}_dedup"

# If input is hg38 liftover to hg19, based on code from topmed_imputation repo
if [ "$build" == "hg38" ]; then
	echo "Input build is hg38. Lifting over to hg19."

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
else
	echo "Input build is already hg19."
	plink \
		--bfile "${out_dir}/tmp_${plink_output}_dedup" \
		--keep-allele-order \
		--make-bed --out "${out_dir}/${plink_output}_dedup_hg19"
fi

# Clean up
rm ${out_dir}/tmp_*
