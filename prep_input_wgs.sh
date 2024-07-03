#!/bin/bash

# Script preps input for local ancestry pipeline, currently
# written to start from a WGS input split by chromosome in
# hg38.

# Get args
chr=17
vcf_file="background/from_dayam_server/topmed_freeze8_phased_sarp/chr${chr}.vcf.gz"
out_dir="data/rfmix_wgs/input"
min_maf=0.03
n_cpus=10

# Remove ID, INFO, FORMAT fields from VCF file
bcftools  annotate -x ID,INFO,FORMAT \
    "$vcf_file" -Oz -o "${out_dir}/tmp_chr${chr}_noinfo.vcf.gz" \
    --threads $nr_threads

# Filter by MAF
bcftools view -i "MAF>0.03" \
    "${out_dir}/chr${chr}_noinfo.vcf.gz" \
    -Oz -o  "${out_dir}/chr${chr}_noinfo_maf${min_maf}.vcf.gz" \
    --threads $nr_threads

# Update variant IDs
snp_id_vcf_file="${out_dir}/chr${chr}_noinfo_maf${min_maf}_snp_id.vcf.gz"
bcftools annotate --set-id '%CHROM:%POS:%REF:%FIRST_ALT' \
    "${out_dir}/chr${chr}_noinfo_maf${min_maf}.vcf.gz" | \
    bgzip -c > $snp_id_vcf_file
tabix $snp_id_vcf_file

# Clean up
rm ${out_dir}/tmp_*
