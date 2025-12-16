#!/bin/bash

# Script requires batched output run with RFMix -co option

lai_dir=$1
out_dir=$2
local_anc_dir=$3
plink_file_name=$4

# Get current code dir
code_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Create PLINK transposed files & local_ancestry_annot_chr#.txt files
for ((chr=1; chr<=22; chr++)); do
    echo "Processing chromosome" $chr
    cat ${lai_dir}/chr${chr}/chr${chr}_local_ancestry.0.Viterbi.txt >> tmp_anc.txt
    n=`wc -l ${lai_dir}/chr${chr}/chr${chr}_local_ancestry_batch1.0.SNPsPerWindow.txt |  tr -s ' ' | cut -f2 -d' '`
    cat ${code_dir}/calc_segment_pos.R | R --vanilla --args $chr $lai_dir $local_anc_dir
    cut -f8 ${lai_dir}/local_ancestry_annot_chr${chr}.txt | sed -e '1d' >> tmp_pos.txt
    for ((i=1; i<=$n; i++)); do
        echo $chr >>  tmp_chr.txt
        echo "0" >> tmp_0.txt
        echo $chr:$i >> tmp_snp.txt
    done
done
paste tmp_chr.txt tmp_snp.txt tmp_0.txt tmp_pos.txt tmp_anc.txt > ${out_dir}/anc.tped
cp ${plink_file_name}.fam ${out_dir}/anc.tfam

#Convert transposed files to a normal format
plink --tfile ${out_dir}/anc --recode 12 --out ${out_dir}/local_ancestry
plink --file ${out_dir}/local_ancestry --make-bed --out ${out_dir}/local_ancestry_final

# Cleanup
rm tmp_*