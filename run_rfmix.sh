#!/bin/bash

#Set parameters
shapeit_haps_file=$1
shapeit_samples_file=$2
rfmix_input_dir=$3

#Extract chromosome number from shapeit file name
chr=`basename $shapeit_haps_file | sed 's/chr//' | sed 's/.haps//'`

#Get list of hapmap SNPs that are in genetic map and also their ref alleles
#Output file: tmp_chr${chr}_tgp_mapped.txt
cut -f2 -d' ' ${rfmix_input_dir}/genetic_map_tgp/chr${chr}.txt | sort -k1 > tmp_chr${chr}_map_sorted.txt
cat ${rfmix_input_dir}/tgp/chr${chr}.impute.legend | grep ^rs | sort -k2 > tmp_chr${chr}_tgp_sorted.txt
cat ${rfmix_input_dir}/tgp/chr${chr}.impute.legend | grep ^rs | sort -k2 | cut -f2 -d' ' | uniq -u > tmp_chr${chr}_tgp_unique.txt
join tmp_chr${chr}_map_sorted.txt tmp_chr${chr}_tgp_unique.txt > tmp_chr${chr}_unique.txt
join -1 1 -2 2 tmp_chr${chr}_unique.txt tmp_chr${chr}_tgp_sorted.txt > tmp_chr${chr}_tgp_mapped.txt

#Get list of admixed SNPs that are in the tmp_map and
#-is not AT/CG
#-do not have allele mismatches
#-determine if ref SNP has changed
#Output file: tmp_chr${chr}_snps_keep.txt
cat $shapeit_haps_file | cut -f2-5 -d ' ' > tmp_chr${chr}_admixed_snps.txt
python3 /home/rfmix_file_creation_scripts/get_admixed_snps.py tmp_chr${chr}_tgp_mapped.txt tmp_chr${chr}_admixed_snps.txt tmp_chr${chr}_snps_keep.txt

#Rewrite the reference population file in the right format and to contain only the SNPs to keep
#Output file: tmp_chr${chr}_ref_haps.txt
sed -e '1d'  ${rfmix_input_dir}/tgp/chr${chr}.impute.legend | cut -f2 -d' '  | sed "s/^/${chr}:/" >  tmp_chr${chr}_ref_snps.txt
paste tmp_chr${chr}_ref_snps.txt  ${rfmix_input_dir}/tgp/chr${chr}.impute.hap > tmp_chr${chr}_ref_select.txt
python3 /home/rfmix_file_creation_scripts/get_ref_haps.py tmp_chr${chr}_snps_keep.txt tmp_chr${chr}_ref_select.txt tmp_chr${chr}_ref_haps.txt

#Rewrite the admixed file in the right format and to contain only the SNPs to keep, and change 0/1 codings where necessary for different ref/alt alleles
#Output file: tmp_chr${chr}_admixed_haps.txt 
python3 /home/rfmix_file_creation_scripts/get_admixed_haps.py tmp_chr${chr}_snps_keep.txt $shapeit_haps_file tmp_chr${chr}_admixed_haps.txt 

#Create the alleles.txt file
#Output file:  chr${chr}_alleles.txt 
paste -d' ' tmp_chr${chr}_admixed_haps.txt tmp_chr${chr}_ref_haps.txt  > tmp_chr${chr}_alleles.txt
sed 's/ //g' tmp_chr${chr}_alleles.txt >  chr${chr}_alleles.txt 

#Create the classes.txt file
#Output file: chr${chr}_classes.txt 
nr_samples=`wc -l $shapeit_samples_file | xargs | cut -f1 -d' '`
let "nr_samples=nr_samples-2"
python3 /home/rfmix_file_creation_scripts/create_classes.py chr${chr}_classes.txt $nr_samples

#Create the snp_locations.txt file
#Output file: chr${chr}_snp_locations.txt
cut -f1 -d' ' tmp_chr${chr}_snps_keep.txt > tmp_chr${chr}_final_pos.txt
cut -f2-3 -d' ' ${rfmix_input_dir}/genetic_map_tgp/chr${chr}.txt | uniq > tmp_chr${chr}_pos_cm.txt
join -1 1 -2 1 tmp_chr${chr}_final_pos.txt tmp_chr${chr}_pos_cm.txt | cut -f2 -d' ' > chr${chr}_snp_locations.txt

#Save the SNP details
#Output file: chr${chr}_local_ancestry_snps.txt
cp tmp_chr${chr}_snps_keep.txt chr${chr}_local_ancestry_snps.txt

#Save the sample information
#Output file: chr${chr}_local_ancestry_samples.txt
cat $shapeit_samples_file | sed -e '1,2d' | cut -f1,2 -d' ' > chr${chr}_local_ancestry_samples.txt

#Run RFMix
RFMix_PopPhased -a chr${chr}_alleles.txt \
                -p chr${chr}_classes.txt \
                -m chr${chr}_snp_locations.txt \
                -co 1 -o chr${chr}_local_ancestry
