#!/bin/bash

chr=$1
begin_pos=$2
end_pos=$3
pheno=$4
data_dir=$5
out_dir=$6

mkdir "$out_dir"

# Get current code dir
code_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Reformat SNP list from RFMix
cut -f1,3,4 -d' ' "${data_dir}/chr${chr}_local_ancestry_snps.txt" > "$data_dir"/snp_info.txt

# Get just IIDs from list of FIDs & IIDs output by RFMix
cut -f2 -d' ' "${data_dir}/chr${chr}_local_ancestry_samples.txt" > "$data_dir"/sample_ids.txt

Rscript ${code_dir}/get_coord.R \
   chr${chr} $begin_pos $end_pos $data_dir $out_dir

new_begin_pos=`cat "$out_dir"/tmp_chr${chr}_begin.txt`
new_end_pos=`cat "$out_dir"/tmp_chr${chr}_end.txt`

python ${code_dir}/create_dose_frames.py \
   "${data_dir}/chr${chr}_local_ancestry.allelesRephased0.txt" \
   "${data_dir}/chr${chr}_local_ancestry.0.Viterbi.txt" \
   "${data_dir}/snp_info.txt" \
   "${data_dir}/sample_ids.txt" \
   $new_begin_pos $new_end_pos \
   "${out_dir}/${pheno}_chr_${chr}_"

# Create second output dir for hg38 liftover
out_dir_hg38="${out_dir//hg19/hg38}"
mkdir $out_dir_hg38

# Do liftover
bash ${code_dir}/liftover_to_hg38.sh $chr \
   "$out_dir"/${pheno}_chr_${chr}_allele_dose.txt  \
   "$out_dir_hg38"

# Update all three dose frames with lifted over coordinates
Rscript ${code_dir}/update_coord.R \
   "$out_dir_hg38"/tmp_chr${chr}_out.bed \
   "$out_dir"/${pheno}_chr_${chr}_ \
   "$out_dir_hg38"/${pheno}_chr_${chr}_

# Cleanup
rm ${out_dir}/tmp_*
rm ${out_dir_hg38}/tmp_*