#!/bin/bash

#Set parameters
if [ "$#" -eq  "0" ]
then
    echo "Usage: ${0##*/} <plink_file_name> <rfmix_results_dir>"
    echo "      <rfmix_input_dir> <admix_input_dir>"
    exit
fi

plink_file_name=$1
rfmix_results_dir=$2
rfmix_input_dir=$3
admix_input_dir=$4

# Get an LD pruned list of SNPs to use from admixed populations
Rscript calc_gwide_ancestry/get_initial_admix_snps.R \
    ${plink_file_name}.bim ${rfmix_results_dir} ${admix_input_dir}

plink --bfile $plink_file_name \
      --keep-allele-order \
      --extract  ${admix_input_dir}/initial_admix_snps.txt \
      --make-bed --out ${admix_input_dir}/initial_admixed

plink --bfile ${admix_input_dir}/initial_admixed \
      --keep-allele-order \
      --geno 0.01 \
      --make-bed --out ${admix_input_dir}/geno_filtered_admixed
plink --bfile ${admix_input_dir}/geno_filtered_admixed \
      --keep-allele-order \
      --indep-pairwise 50 10 0.1 \
      --out ${admix_input_dir}/ld
plink --bfile ${admix_input_dir}/geno_filtered_admixed \
      --keep-allele-order \
      --extract ${admix_input_dir}/ld.prune.in \
      --make-bed --out ${admix_input_dir}/ld_filtered_admixed

# Create tped file from reference panels and convert it to a BED file
rm "${rfmix_input_dir}/working/ref.tped"

for ((chr=1; chr<=22; chr++));
do
    sed -e '1d' ${rfmix_input_dir}/tgp/chr${chr}.impute.legend > \
        ${rfmix_input_dir}/working/chr${chr}_snp_info.txt

    paste ${rfmix_input_dir}/working/chr${chr}_snp_info.txt \
          ${rfmix_input_dir}/tgp/chr${chr}.impute.hap > \
          ${rfmix_input_dir}/working/chr${chr}.hap

    grep "^$chr\t" ${admix_input_dir}/ld_filtered_admixed.bim | \
        cut -f4 > ${admix_input_dir}/chr${chr}.keep

    python calc_gwide_ancestry/create_ref_tped_files.py \
        $chr $admix_input_dir ${rfmix_input_dir}/working

    cat ${rfmix_input_dir}/working/chr${chr}.tped >> ${rfmix_input_dir}/working/ref.tped
done

# Make TGP ref .fam file and then PLINK file
n=`wc -l ${rfmix_input_dir}/tgp/chr22.impute.hap.indv | xargs | cut -f1 -d' '`
rm ${admix_input_dir}/tmp_c.txt

for ((i=1; i<=$n; i++)); do
    echo "0" >> ${admix_input_dir}/tmp_c.txt
done

paste  ${rfmix_input_dir}/tgp/chr22.impute.hap.indv \
       ${rfmix_input_dir}/tgp/chr22.impute.hap.indv \
       ${admix_input_dir}/tmp_c.txt \
       ${admix_input_dir}/tmp_c.txt \
       ${admix_input_dir}/tmp_c.txt \
       ${admix_input_dir}/tmp_c.txt > ${rfmix_input_dir}/working/ref.tfam
plink --tfile ${rfmix_input_dir}/working/ref \
      --keep-allele-order \
      --make-bed --out ${rfmix_input_dir}/working/b_ref

# Before merge, use PLINK2 to set variant ids to chr:pos and flip A2
# alleles in dataset to match reference. Many of the remaining warnings
# at this step should be fixed by flipping after first merge attempt
plink2 --bfile ${rfmix_input_dir}/working/b_ref \
    --set-all-var-ids @:# --make-bed \
    --out ${rfmix_input_dir}/working/b_ref_chrpos
plink --bfile ${admix_input_dir}/ld_filtered_admixed \
    --keep-allele-order \
    --a2-allele ${rfmix_input_dir}/working/b_ref_chrpos.bim 6 2 \
    --make-bed --out ${admix_input_dir}/ld_filtered_admixed_a2_ref

# Merge the reference and admixed files
plink --bfile ${rfmix_input_dir}/working/b_ref_chrpos --allow-no-sex \
      --keep-allele-order \
      --bmerge ${admix_input_dir}/ld_filtered_admixed_a2_ref \
      --make-bed --out ${rfmix_input_dir}/working/dummy_merge

if [ -e "${rfmix_input_dir}/working/dummy_merge-merge.missnp" ]
then
    mv ${rfmix_input_dir}/working/dummy_merge-merge.missnp ${rfmix_input_dir}/working/flip_snps.txt
    plink --bfile ${admix_input_dir}/ld_filtered_admixed_a2_ref \
          --keep-allele-order \
          --flip ${rfmix_input_dir}/working/flip_snps.txt \
          --make-bed --out ${rfmix_input_dir}/working/flip_admixed
else
    mv ${rfmix_input_dir}/working/dummy_merge.bed ${rfmix_input_dir}/working/flip_admixed.bed
    mv ${rfmix_input_dir}/working/dummy_merge.bim ${rfmix_input_dir}/working/flip_admixed.bim
    mv ${rfmix_input_dir}/working/dummy_merge.fam ${rfmix_input_dir}/working/flip_admixed.fam
fi
plink --bfile ${rfmix_input_dir}/working/b_ref_chrpos \
      --keep-allele-order \
      --bmerge ${rfmix_input_dir}/working/flip_admixed \
      --allow-no-sex \
      --make-bed --out ${rfmix_input_dir}/working/dummy_merge
if [ -e "${rfmix_input_dir}/working/dummy_merge-merge.missnp" ]
then
    plink --bfile  ${rfmix_input_dir}/working/flip_admixed \
          --keep-allele-order \
          --remove ${rfmix_input_dir}/working/dummy_merge-merge.missnp \
          --make-bed --out ${rfmix_input_dir}/working/flip_fixed_admixed
    plink --bfile ${rfmix_input_dir}/working/b_ref_chrpos \
      --keep-allele-order \
      --bmerge ${rfmix_input_dir}/working/flip_fixed_admixed \
      --allow-no-sex \
      --make-bed --out ${rfmix_input_dir}/working/dummy_merge
    rm ${rfmix_input_dir}/working/dummy_merge-merge.missnp
fi

# Move to admixutre input folder
mv ${rfmix_input_dir}/working/dummy_merge.bed ${admix_input_dir}/merged.bed
mv ${rfmix_input_dir}/working/dummy_merge.bim ${admix_input_dir}/merged.bim
mv ${rfmix_input_dir}/working/dummy_merge.fam ${admix_input_dir}/merged.fam
mv ${rfmix_input_dir}/working/dummy_merge.log ${admix_input_dir}/merged.log
