#!/bin/bash

# Outputs dose frames in hg19 and hg38 based on admixture mapping peak
# regions written out by summary report. Regions are contiguous segments
# with P < 1E-3 around significant peaks.

out_dir_prefix="$1"
rfmix_results_dir="$2"
admix_dir="$3"

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
        $chr $start $end $pheno "${rfmix_results_dir}/chr${chr}" $out_dir
done
