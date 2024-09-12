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

# Set VCF output
vcf_output=`basename $vcf_input | sed 's/.vcf.gz//'`

# Remove ID, INFO, FORMAT fields from VCF file
bcftools  annotate -x ID,^INFO/R2,INFO/MAF,FORMAT \
     "$vcf_input" -Oz -o "${out_dir}/tmp_${vcf_output}_noinfo.vcf.gz" \
     --threads $nr_threads

# Because bcftools view with MAF filter is silent if MAF tag absent,
# use bcftools to fill MAF values
bcftools +fill-tags "${out_dir}/tmp_${vcf_output}_noinfo.vcf.gz" \
     -Oz -o "${out_dir}/tmp_${vcf_output}_noinfo_maf.vcf.gz"

# Filter by MAF
bcftools view -i "MAF>${min_maf}" \
     "${out_dir}/tmp_${vcf_output}_noinfo_maf.vcf.gz" \
     -Oz -o  "${out_dir}/tmp_${vcf_output}_noinfo_maf${min_maf}.vcf.gz" \
     --threads $nr_threads

# Update variant IDs
snp_id_vcf_file="${out_dir}/tmp_${vcf_output}_noinfo_maf${min_maf}_snp_id.vcf.gz"
bcftools annotate --set-id '%CHROM:%POS:%REF:%FIRST_ALT' \
     "${out_dir}/tmp_${vcf_output}_noinfo_maf${min_maf}.vcf.gz" | \
     bgzip -c > $snp_id_vcf_file
tabix $snp_id_vcf_file

# Filter to just samples of interest
snp_id_vcf_file_filt="${out_dir}/${vcf_output}_noinfo_maf${min_maf}_snp_id_filt.vcf.gz"
bcftools view -S "$samp_ids" \
    $snp_id_vcf_file \
    -Oz -o $snp_id_vcf_file_filt
tabix $snp_id_vcf_file_filt

# Clean up
rm ${out_dir}/tmp_*
