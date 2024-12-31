#!/bin/bash

# TO NOTE: update the paths below to set the chromosome (chr), windows
# (windows, list of "start stop" in hg19), phenotype for analysis (pheno)
# the input data directory (data_dir), and the output directory prefix
# (out_dir_prefix). Outputs dose frames in hg19 and hg38. Input admixture
# mapping peaks are in hg19, files written out by summary report.

# TO NOTE: activate crossmap-bcftools-osx64 before running.

# Scripts assume RFMix output file names are as output by RFMix (needs to
# be run with -o option):
    # chr#_local_ancestry.0.Viterbi.txt
    # chr#_local_ancestry.allelesRephased0.txt
    # chr#_local_ancestry_snps.txt
    # chr#_local_ancestry_samples.txt


# Set inputs
proj_dir="${RKJCOLLAB}/Collabs/ortega"
out_dir_prefix="${proj_dir}/data/laaa_wgs_new_2/dose_frames"
data_dir="${proj_dir}/data/rfmix_wgs_new_2/output_o"
# pheno=maxFEV1_FVC  # switch pheno to run for each region
# chr=10  # switch chr to run for each region
pheno=maxFVC
chr=5
# chr=8
# chr=17

# Region file to be filled in with above phenotype and chromosome, made by admixture
# mapping summary report with one peak and its contiguous region
region="${proj_dir}/data/pheno_analysis/output/${pheno}_admixture_peak_contig_region_chr${chr}.txt"

# Get region start/stop in hg19
start=$(awk 'NR==2 {print $3}' $region | cut -d'-' -f1)
end=$(awk 'END {print $3}' $region | cut -d'-' -f2)

# Create dose frames for each contiguous region
mkdir "$out_dir_prefix"
out_dir="${out_dir_prefix}/region_hg19_${pheno}_chr${chr}"

bash make_dose_frames/create_dose_frames.sh \
    $chr $start $end $pheno "${data_dir}/chr${chr}" $out_dir
