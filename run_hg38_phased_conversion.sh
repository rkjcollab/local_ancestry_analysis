#!/bin/bash

# vcf_file=$1
# min_maf=$2
# min_rsq=$3

# TODO: need to activate conda env before running
# TO NOTE: SDS modifying to accept WGS file instead of imputed file
# TODO: write now, not keeping any INFO fields, but if was going to
# use for imputed files, would need to edit to keep R2 at step of MAF
# TO NOTE: run from top level of local_anc_afr_eur repo

# TODO: vcf file should be after filtering, with SNP IDs
vcf_file=$1
out_dir=$2
n_cpus=$3
code_dir=$4

file_id=`basename $vcf_file | sed 's/.vcf.gz//'`

# #Apply the Rsq filter (if zero all the SNPs will be retained)
# TODO: skipping for now, revisit to edit to work for imputed files
rsq_vcf_file=$vcf_file

# #Get lifted over coordinates of variants 
zcat $rsq_vcf_file | grep -v ^# | cut -f1 > "${out_dir}/tmp_${file_id}_chr_col.txt"
zcat $rsq_vcf_file | grep -v ^# | cut -f2 > "${out_dir}/tmp_${file_id}_pos_col.txt"
zcat $rsq_vcf_file | grep -v ^# | cut -f3 > "${out_dir}/tmp_${file_id}_snp_col.txt"
chr=`head -1 ${out_dir}/tmp_${file_id}_chr_col.txt | sed 's/chr//'`
paste "${out_dir}/tmp_${file_id}_chr_col.txt" \
   "${out_dir}/tmp_${file_id}_pos_col.txt tmp_${file_id}_pos_col.txt" \
   "${out_dir}/tmp_${file_id}_snp_col.txt" > "${out_dir}/tmp_chr${chr}_in.bed"

# Make a BED file: chr, pos (0-based), end pos (1-based), id
# TODO: alt option if the above is slow
# bcftools query -f'%CHROM\t%POS0\t%END\t%ID\n' file.bcf

CrossMap.py bed "${code_dir}/shapeit_formatting_scripts/hg38ToHg19.over.chain" \
            "${out_dir}/tmp_chr${chr}_in.bed"  \
            "${out_dir}/tmp_chr${chr}_out.bed"

#Create bed file to annotate the hg19 position
cat "${code_dir}/shapeit_formatting_scripts/create_hg19_annot_bed.R" | R --vanilla --args $chr \
   "${out_dir}/tmp_chr${chr}_out.bed" \
   "${out_dir}/tmp_chr${chr}_hg19_annot.txt"
bgzip "${out_dir}/tmp_chr${chr}_hg19_annot.txt"
tabix -s1 -b2 -e3 "${out_dir}/tmp_chr${chr}_hg19_annot.txt.gz"

# #Annotate the hg19 column
# hg19_annot_vcf_file="${out_dir}/tmp_${file_id}_min_maf_with_snp_id_min_rsq_hg19_annot.vcf"
# bcftools annotate \
#   -a  tmp_chr${chr}_hg19_annot.txt.gz \
#   -c CHROM,FROM,TO,REF,ALT,HG19 \
#   -h <(echo '##INFO=<ID=HG19,Number=1,Type=Integer,Description="hg19 position">') \
#   $rsq_vcf_file > $hg19_annot_vcf_file

# #Filter out variants without an hg19 annotation
# vcftools --gzvcf $hg19_annot_vcf_file \
#    --get-INFO HG19 \
#    --stdout > tmp_${file_id}_hg19.txt
# grep -v "\?$" tmp_${file_id}_hg19.txt | \
#    cut -f1-4 | \
#    sed -e $'s/\t/:/g' | sed -e '1d' > \
#    tmp_${file_id}_hg19_variants.txt
# hg19_filtered_vcf_file=tmp_${file_id}_hg19_only.vcf
# vcftools --gzvcf $hg19_annot_vcf_file \
#    --snps tmp_${file_id}_hg19_variants.txt \
#    --recode --recode-INFO HG19 \
#    --stdout > $hg19_filtered_vcf_file

# #Update POS column with hg19 position
# hg19_pos_vcf_file=tmp_${file_id}_hg19_pos.vcf
# grep "^#" $hg19_filtered_vcf_file > $hg19_pos_vcf_file
# grep -v "^#" ${hg19_filtered_vcf_file} | cut -f1 > ${hg19_filtered_vcf_file}.c1 #chr
# vcftools --vcf $hg19_filtered_vcf_file \
#    --get-INFO HG19 \
#    --stdout | cut -f5 | sed -e '1d' > ${hg19_filtered_vcf_file}.c2
# grep -v "^#" $hg19_filtered_vcf_file | cut -f 3- >  ${hg19_filtered_vcf_file}.c3_onwards
# paste ${hg19_filtered_vcf_file}.c1 \
#    ${hg19_filtered_vcf_file}.c2 \
#    ${hg19_filtered_vcf_file}.c3_onwards >> $hg19_pos_vcf_file

# #Sort by position
# hg19_sorted_vcf_file=tmp_${file_id}_hg19_sorted.vcf
# cat $hg19_pos_vcf_file | vcf-sort > $hg19_sorted_vcf_file

# #Convert the hg19 VCF to a shapeit format haps file chr${chr}.haps
# grep ^# -v $hg19_sorted_vcf_file > ${hg19_sorted_vcf_file}.genos
# cut -f1 ${hg19_sorted_vcf_file}.genos > ${hg19_sorted_vcf_file}.haps.c1 #chr
# cut -f3 ${hg19_sorted_vcf_file}.genos > ${hg19_sorted_vcf_file}.haps.c2 #snp_id
# cut -f2 ${hg19_sorted_vcf_file}.genos > ${hg19_sorted_vcf_file}.haps.c3 #pos
# cut -f4 ${hg19_sorted_vcf_file}.genos > ${hg19_sorted_vcf_file}.haps.c4 #ref
# cut -f5 ${hg19_sorted_vcf_file}.genos > ${hg19_sorted_vcf_file}.haps.c5 #alt
# #Sometimes more than just GT is encoded in a VCF file so extract only the genotypes
# vcftools --vcf $hg19_sorted_vcf_file \
#    --extract-FORMAT-info GT \
#    --stdout | \
#    cut -f 3- | sed -e 's/|/ /g' | \
#    expand -t1 | \
#    sed -e '1d' > \
#    ${hg19_sorted_vcf_file}.haps.c6_onwards #genotypes
# paste -d' ' ${hg19_sorted_vcf_file}.haps.c1 \
#     ${hg19_sorted_vcf_file}.haps.c2 \
#     ${hg19_sorted_vcf_file}.haps.c3 \
#     ${hg19_sorted_vcf_file}.haps.c4 \
#     ${hg19_sorted_vcf_file}.haps.c5 \
#     ${hg19_sorted_vcf_file}.haps.c6_onwards > chr${chr}.haps

# #Create the samples file chr${chr}.samples
# sample_line=`grep \#CHROM $hg19_sorted_vcf_file | cut -f 10-` 
# echo "ID2" > chr${chr}.samples
# echo "0" >> chr${chr}.samples
# for value in $sample_line
# do
#     echo $value >> chr${chr}.samples
# done

