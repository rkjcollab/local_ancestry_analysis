#!/bin/bash

# TO NOTE: update the paths below to set the path to the PLINK1.9 file set used
# as input to RFMIX with all variant IDs set to chr:pos:ref:alt (plink_file_name),
# the path to RFMix results split by chromosome (rfmix_results dir), path to
# the unzipped RFMix input directory containing reference files (rfmix_input_dir),
# the path to the location to save ADMIXTURE inputs and outputs (admix_dir), and
# path where want to save combined genome wide ancestry results (out_dir).

# Scripts assume RFMix output file names are as output by RFMix:
    # chr#_local_ancestry.0.Viterbi.txt
    # chr#_local_ancestry.allelesRephased0.txt
    # chr#_local_ancestry_snps.txt
    # chr#_local_ancestry_samples.txt
    # chr#_local_ancestry_batch1.0.SNPsPerWindow.txt (if found in dir, means
        # RFMix was run with -co option)

plink_file_name="background/from_dayam_server/Victor_Ortega_analyses_MPB/SARP_CSGA_merge/SARP123_CSGA_merged_final"
rfmix_results_dir="data/rfmix_unimp/output_co"
rfmix_input_dir="/Users/slacksa/temp_ortega/rfmix_input"  # genetic map & other unchanging inputs
admix_dir="data/admixture_unimp"
out_dir="data/gwide_anc_unimp"

# Calculate RFMIX genome-wide ancestry
    # Third optional argument is SNPs per window file. Only need if RFMix
    # was not run in batches as that means code won't auto find the file.
Rscript code/calc_gwide_ancestry/calc_rfmix_gwide_ancestry.R \
    ${plink_file_name}.fam $rfmix_results_dir
    # ${plink_file_name}.fam $rfmix_results_dir "local_ancestry.0.SNPsPerWindow.txt"

# Calculate ADMIXTURE genome-wide ancestry
# Create working dir below rfmix_input_dir, delete at end
rm -r "${rfmix_input_dir}/working"
mkdir "${rfmix_input_dir}/working"
# Create admixture dirs
mkdir "${admix_dir}"
admix_input_dir="${admix_dir}/input"
mkdir "${admix_input_dir}"
admix_out_dir="${admix_dir}/output"
mkdir "${admix_out_dir}"

# Make admixture inputs
bash code/calc_gwide_ancestry/create_admixture_input.sh \
    $plink_file_name $rfmix_results_dir $rfmix_input_dir $admix_input_dir

# Run admixture with k = 2
admixture "${admix_input_dir}/merged.bed" 2

# Move results from working dir to output dir
mv merged.* $admix_out_dir
nr_samples=`wc -l ${plink_file_name}.fam | xargs | cut -f1 -d' '`
Rscript code/calc_gwide_ancestry/calc_admixture_gwide_ancestry.R \
    "${plink_file_name}.fam" "${admix_input_dir}/merged.fam" \
    $admix_input_dir $admix_out_dir \
    "${rfmix_input_dir}/igsr-1000_genomes_phase_3_release.tsv"

#Create a combined genome-wide output file
mkdir "${out_dir}"
Rscript code/calc_gwide_ancestry/merge_gwide_ancestry.R \
    $rfmix_results_dir $admix_out_dir $out_dir

# Clean up
# rm -r "${rfmix_input_dir}/working"
