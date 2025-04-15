#!/bin/bash

if [ "$#" -eq  "0" ]
then
   echo "Usage: ${0##*/} <path_to_vcf_input> <output_dir_prefix> <min_maf>"
   echo "      <nr_threads> <sample_id_list>"
   echo "Script preps WGS or imputed input for local ancestry pipeline. Input"
   echo "should be a single-chromosome VCF file in either hg19 or hg38."
   exit
fi

# Get args
vcf_input=$1
out_dir=$2
min_maf=$3
nr_threads=$4
samp_ids=$5

vcf_prefix=$(basename "$vcf_input")
vcf_prefix="${vcf_prefix%.gz}"  # removes .gz if present
vcf_prefix="${vcf_prefix%.vcf}"  # removes .vcf

# Remove ID, INFO, FORMAT fields from VCF file
     # including format removes all tags except for GT
bcftools  annotate -x ID,^INFO/R2,INFO/MAF,FORMAT \
     "$vcf_input" -Oz -o "${out_dir}/${vcf_prefix}_noinfo.vcf.gz" \
     --threads $nr_threads

# Filter to just samples of interest
bcftools view -S "$samp_ids" \
    "${out_dir}/${vcf_prefix}_noinfo.vcf.gz" \
    -Oz -o "${out_dir}/tmp_${vcf_prefix}_noinfo_filt.vcf.gz"

# Update MAF calculations after sample filtering
bcftools +fill-tags "${out_dir}/tmp_${vcf_prefix}_noinfo_filt.vcf.gz" \
     -Oz -o "${out_dir}/tmp_${vcf_prefix}_noinfo_filt_tags.vcf.gz" \
     --threads $nr_threads -- -t MAF

# Filter by MAF
bcftools view -i "MAF>${min_maf}" "${out_dir}/tmp_${vcf_prefix}_noinfo_filt_tags.vcf.gz" \
     -Oz -o  "${out_dir}/tmp_${vcf_prefix}_noinfo_filt_tags_maf${min_maf}.vcf.gz" \
     --threads $nr_threads

# Update variant IDs
snp_id_vcf_file="${out_dir}/${vcf_prefix}_noinfo_filt_tags_maf${min_maf}_snp_id"
bcftools annotate --set-id '%CHROM:%POS:%REF:%FIRST_ALT' \
    "${out_dir}/tmp_${vcf_prefix}_noinfo_filt_tags_maf${min_maf}.vcf.gz" \
     --threads $nr_threads | \
     bgzip -c > "${snp_id_vcf_file}.vcf.gz"
tabix "${snp_id_vcf_file}.vcf.gz"

# Clean up
rm ${out_dir}/tmp_*
