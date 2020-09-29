#!/bin/bash

vcf_file=$1
min_maf=$2
min_rsq=$3
file_id=`basename $vcf_file | sed 's/.vcf.gz//'`

#First create a potentially smaller file by applying the MAF filter
if [ "$min_maf" == "0" ] 
then
   maf_vcf_file=$vcf_file
else 
   maf_vcf_file=tmp_${file_id}_min_maf${min_maf}.vcf.gz
   vcftools --gzvcf $vcf_file \
      --maf $min_maf --recode --recode-INFO R2 \
      --stdout | bgzip -c > $maf_vcf_file
   tabix $maf_vcf_file
fi

#Update ID field to chr:pos:ref:alt, as this field is often set to '.'
#(at least that is the case in TOPMed files)
#A SNP ID is required in order for the downstream steps to work
#The original hg38 position is also extracted from the SNP ID in later steps
snp_id_vcf_file=tmp_${file_id}_min_maf_with_snp_id.vcf.gz
bcftools annotate --set-id +'%CHROM:%POS:%REF:%FIRST_ALT' $maf_vcf_file | \
   bgzip -c > $snp_id_vcf_file
tabix $snp_id_vcf_file

#Apply the Rsq filter (if zero all the SNPs will be retained)
vcftools --gzvcf $snp_id_vcf_file \
   --get-INFO R2 \
   --stdout > tmp_${file_id}_r2.txt
if [ "$min_rsq" == "0" ] 
then
   touch tmp_${file_id}_r2_failed_variants.txt #create an empty file 
else 
   cat filter_rsq.R | R --vanilla --args tmp_${file_id}_r2.txt \
      tmp_${file_id}_r2_failed_variants.txt \
      $min_rsq
fi
if [ -s "tmp_${file_id}_r2_failed_variants.txt" ]
then
   rsq_vcf_file=tmp_${file_id}_min_maf_with_snp_id_min_rsq.vcf.gz
   vcftools --gzvcf $snp_id_vcf_file \
      --exclude tmp_${file_id}_r2_failed_variants.txt \
      --recode --stdout | bgzip -c > $rsq_vcf_file
else
   rsq_vcf_file=$snp_id_vcf_file
fi

#Get lifted over coordinates of variants 
zcat $rsq_vcf_file | grep -v ^# | cut -f1 > tmp_${file_id}_chr_col.txt
zcat $rsq_vcf_file | grep -v ^# | cut -f2 > tmp_${file_id}_pos_col.txt
zcat $rsq_vcf_file | grep -v ^# | cut -f3 > tmp_${file_id}_snp_col.txt
chr=`head -1 tmp_chr22_chr_col.txt | sed 's/chr//'`
paste tmp_${file_id}_chr_col.txt \
   tmp_${file_id}_pos_col.txt tmp_${file_id}_pos_col.txt \
   tmp_${file_id}_snp_col.txt > tmp_chr${chr}_in.bed 
CrossMap.py bed /home/shapeit_formatting_scripts/hg38ToHg19.over.chain \
            tmp_chr${chr}_in.bed  \
            tmp_chr${chr}_out.bed

#Create bed file to annotate the hg19 position
cat /home/shapeit_formatting_scripts/create_hg19_annot_bed.R | R --vanilla --args $chr \
   tmp_chr${chr}_out.bed \
   tmp_chr${chr}_hg19_annot.txt
bgzip tmp_chr${chr}_hg19_annot.txt
tabix -s1 -b2 -e3 tmp_chr${chr}_hg19_annot.txt.gz

#Annotate the hg19 column
hg19_annot_vcf_file=tmp_${file_id}_min_maf_with_snp_id_min_rsq_hg19_annot.vcf
bcftools annotate \
  -a  tmp_chr${chr}_hg19_annot.txt.gz \
  -c CHROM,FROM,TO,REF,ALT,HG19 \
  -h <(echo '##INFO=<ID=HG19,Number=1,Type=Integer,Description="hg19 position">') \
  $rsq_vcf_file > $hg19_annot_vcf_file

#Filter out variants without an hg19 annotation
#Also update POS column with hg19 position
hg19_filtered_vcf_file=tmp_${file_id}_hg19_only.vcf
grep ^# $hg19_annot_vcf_file > $hg19_filtered_vcf_file #Add the header lines
grep "HG19=" $hg19_annot_vcf_file > ${hg19_filtered_vcf_file}.genos #Keep only entries with HG19 positions
cut -f1 ${hg19_filtered_vcf_file}.genos > ${hg19_filtered_vcf_file}.c1 #chr
cut -f8 ${hg19_filtered_vcf_file}.genos | sed 's/HG19=//' > ${hg19_filtered_vcf_file}.c2 #hg19 position
cut -f 3- ${hg19_filtered_vcf_file}.genos >  ${hg19_filtered_vcf_file}.c3_onwards
paste ${hg19_filtered_vcf_file}.c1 \
   ${hg19_filtered_vcf_file}.c2 \
   ${hg19_filtered_vcf_file}.c3_onwards >> $hg19_filtered_vcf_file

#Sort by position
hg19_sorted_vcf_file=tmp_${file_id}_hg19_sorted.vcf
cat $hg19_filtered_vcf_file | vcf-sort > $hg19_sorted_vcf_file

#Convert the hg19 VCF to a shapeit format haps file chr${chr}.haps
grep ^# -v $hg19_sorted_vcf_file > tmp_${hg19_sorted_vcf_file}.genos
cut -f1 tmp_${hg19_sorted_vcf_file}.genos > tmp_${hg19_sorted_vcf_file}.haps.c1 #chr
cut -f3 tmp_${hg19_sorted_vcf_file}.genos > tmp_${hg19_sorted_vcf_file}.haps.c2 #snp_id
cut -f2 tmp_${hg19_sorted_vcf_file}.genos > tmp_${hg19_sorted_vcf_file}.haps.c3 #pos
cut -f4 tmp_${hg19_sorted_vcf_file}.genos > tmp_${hg19_sorted_vcf_file}.haps.c4 #ref
cut -f5 tmp_${hg19_sorted_vcf_file}.genos > tmp_${hg19_sorted_vcf_file}.haps.c5 #alt
#Sometimes more than just GT is encoded in a VCF file so extract only the genotypes
vcftools --vcf $hg19_sorted_vcf_file \
   --extract-FORMAT-info GT \
   --recode --stdout | \
   cut -f 3- | sed -e 's/|/ /g' | expand -t1 > \ #replace pipe with spaces and tabs with spaces
   tmp_${hg19_sorted_vcf_file}.haps.c6_onwards #genotypes
paste -d' ' tmp_${hg19_sorted_vcf_file}.haps.c1 \
    tmp_${hg19_sorted_vcf_file}.haps.c2 \
    tmp_${hg19_sorted_vcf_file}.haps.c3 \
    tmp_${hg19_sorted_vcf_file}.haps.c4 \
    tmp_${hg19_sorted_vcf_file}.haps.c5 \
    tmp_${hg19_sorted_vcf_file}.haps.c6_onwards > chr${chr}.haps


#Create the samples file chr${chr}.samples
sample_line=`grep \#CHROM $hg19_sorted_vcf_file | cut -f 10-` 
touch chr${chr}.samples
for value in $sample_line
do
    echo $value >> chr${chr}.samples
done

