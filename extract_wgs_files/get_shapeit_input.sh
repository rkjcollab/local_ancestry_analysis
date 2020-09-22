#!/bin/bash

#Set parameters
chr=$1
if [ "$#" -eq  "0" ]
then
    echo "Usage: ${0##*/} <chr_nr>"
    exit
fi

#Get list of samples to extract
if [ ! -f SARP_AA.txt ]; then
    grep NWD /Users/dayam/Remote/OneDrive\ -\ The\ University\ of\ Colorado\ Denver/Ortega/analysis/input_phenotypes_genotypes_fine_mapping/id_lists/AA*WGS*txt | cut -f2 -d ":" | sort | uniq -d > tmp_SARP_AA.txt
    paste tmp_SARP_AA.txt tmp_SARP_AA.txt > SARP_AA.txt
fi

#Create shapeit direcetory
if [ ! -d ./shapeit_input ]; then
    mkdir shapeit_input
fi

#Get just cleaned SARP AA
plink --vcf "/Volumes/Promise Pegasus/topmed_freeze6a_sarp/chr${chr}.vcf.gz" --keep SARP_AA.txt --make-bed --out tmp_chr${chr}
plink --bfile tmp_chr${chr} --hwe 0.000001 --geno 0.01 --maf 0.005 --snps-only --make-bed --out tmp_hwe_miss_filtered_chr${chr}

#Update SNP names to chr:postion so as to remove duplicate SNPs
cat update_map.R | R --vanilla --args tmp_hwe_miss_filtered_chr${chr}.bim $chr
mv tmp_new_chr${chr}.bim tmp_hwe_miss_filtered_chr${chr}.bim

#Remove SNPs with duplicate positions
cat get_dupl_SNPs.R | R --vanilla --args tmp_hwe_miss_filtered_chr${chr}.bim $chr
plink --bfile tmp_hwe_miss_filtered_chr${chr} \
      --exclude tmp_dupl_chr${chr}.txt \
      --make-bed --out tmp_no_dupl_variants_chr${chr}

#Get input file for crossover
cat tmp_no_dupl_variants_chr${chr}.bim | cut -f1 | sed 's/^/chr/' > tmp_c1_chr${chr}.txt
cat tmp_no_dupl_variants_chr${chr}.bim | cut -f4 > tmp_c2_chr${chr}.txt
cat tmp_no_dupl_variants_chr${chr}.bim | cut -f4 > tmp_c3_chr${chr}.txt
cat tmp_no_dupl_variants_chr${chr}.bim | cut -f2 > tmp_c4_chr${chr}.txt
paste  tmp_c1_chr${chr}.txt \
       tmp_c2_chr${chr}.txt \
       tmp_c3_chr${chr}.txt \
       tmp_c4_chr${chr}.txt \
       >  tmp_in_chr${chr}.bed

#Do crossover
CrossMap.py bed hg38ToHg19.over.chain \
            tmp_in_chr${chr}.bed  \
            tmp_out_chr${chr}.bed

#Extract only those SNPs that were successfully cross-overed (on the same chromosome)
grep  "^chr${chr}\t" tmp_out_chr${chr}.bed | cut -f4 > tmp_snp_keep_chr${chr}.txt
plink --bfile tmp_no_dupl_variants_chr${chr} \
      --extract tmp_snp_keep_chr${chr}.txt \
      --make-bed --out tmp_clean_chr${chr}

#Update bim file positions and SNP names
cat liftover_map.R | R --vanilla --args tmp_clean_chr${chr}.bim $chr
mv tmp_lifted_chr${chr}.bim tmp_clean_chr${chr}.bim

#Run PLINK again since order of SNPs might have changed with liftover
#Then you can also identify duplicate SNPs that may have been introduced with liftover
plink --bfile tmp_clean_chr${chr} \
      --make-bed --out tmp_ordered_chr${chr}
cut -f2 tmp_ordered_chr${chr}.bim | uniq -d > tmp_dupl_chr${chr}.txt
plink --bfile tmp_ordered_chr${chr} \
      --exclude tmp_dupl_chr${chr}.txt \
      --make-bed --out ./shapeit_input/chr${chr}

rm tmp_*chr${chr}.*
