#!/bin/bash

# Set inputs
pheno_analysis_dir=$1
input_dir=$2
out_dir=$3
perc_variance_exp="0.995" # set 0.995 to match PMID 34762840, although Michelle's code had 0.998

# Get current code dir
code_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Make temp working dir
work_dir="${pheno_analysis_dir}/working"
mkdir $work_dir

# Calculate effective tests, using input made by script create_plink_input_local
echo -e chr'\t'm_eff'\t'm > "${out_dir}/m_eff_${perc_variance_exp}.txt"
for ((chr=1; chr<=22; chr++)); do
    grep "^$chr\t" "${input_dir}/local_ancestry.map" | cut -f2 > "${work_dir}/snps.txt"
    plink --file "${input_dir}/local_ancestry" --extract "${work_dir}/snps.txt" \
          --recodeA --out "${work_dir}/chr${chr}"
    Rscript ${code_dir}/calc_m_eff.R \
        $chr $work_dir $out_dir $perc_variance_exp
done
m_eff=`cut -f2 "${out_dir}/m_eff_${perc_variance_exp}.txt" | sed 's/ //g' | sed -e '1d' | paste -sd+ - | bc`
echo $m_eff

# Clean up
rm -r $work_dir