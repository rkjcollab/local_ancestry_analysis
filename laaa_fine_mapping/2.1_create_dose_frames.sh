#!/bin/bash

# TO NOTE: update the paths below to set the project directory (proj_dir),
# the output directory (out_dir_prefix), the input data directory (data_dir),
# and the directory containing admixture mapping peak regions in hg19 
# written out by summary report (admix_dir).

# Script outputs dose frames in hg19 and hg38.

# Scripts assume RFMix output file names are as output by RFMix (needs to
# be run with -o option):
    # chr#_local_ancestry.0.Viterbi.txt
    # chr#_local_ancestry.allelesRephased0.txt
    # chr#_local_ancestry_snps.txt
    # chr#_local_ancestry_samples.txt

# Set inputs
# Test WGS
proj_dir="${RKJCOLLAB}/Collabs/ortega"
out_dir_prefix="${proj_dir}/data/pipeline_test_data/laaa_fine_mapping/laaa_wgs/dose_frames"
data_dir="${proj_dir}/data/rfmix_wgs_new_2/output_o"
admix_dir="${proj_dir}/data/pipeline_test_data/laaa_fine_mapping/pheno_analysis/output"

# Loop over peak admixture mapping regions and make dose frames for each
admix_file_list=$(ls $admix_dir | grep -E "region_chr[0-9]+\.txt")
for admix_file in $admix_file_list; do
    echo "Making dose frames for region file ${admix_file}."

    [[ $admix_file =~ (.*)_admixture_peak_contig_region_chr(.*)\.txt ]];
    pheno=${BASH_REMATCH[1]}
    chr=${BASH_REMATCH[2]}

    # Get region start/stop in hg19
    start=$(awk 'NR==2 {print $3}' "${admix_dir}/${admix_file}" | cut -d'-' -f1)
    end=$(awk 'END {print $3}' "${admix_dir}/${admix_file}" | cut -d'-' -f2)

    # Create dose frames for each contiguous region
    mkdir "$out_dir_prefix"
    out_dir="${out_dir_prefix}/region_hg19_${pheno}_chr${chr}"

    bash make_dose_frames/create_dose_frames.sh \
        $chr $start $end $pheno "${data_dir}/chr${chr}" $out_dir
done
