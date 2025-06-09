
#!/bin/bash

# This is currently called inside of the limactl amd64.
# Set paths and options in config.yml file. Paths are automatically 
# relative to snakefile location.

# Uncomment step want to run through
# step="submit_initial_input"  # pre-imp QC
step="submit_fix_strands"  # submit imp

# Set RKJCOLLAB since not present in lima
RKJCOLLAB="/Users/slacksa/Library/CloudStorage/OneDrive-TheUniversityofColoradoDenver"

# TOPMed
echo "TOPMed"
apptainer exec \
    --writable-tmpfs \
    --bind /Users/slacksa/repos/imputation_snakemake:/repo \
    --bind /Users/slacksa/repos/local_ancestry_analysis/imputation:/proj_repo \
    --bind ${RKJCOLLAB}/Collabs/ortega/background/from_dayam_server/Victor_Ortega_analyses_MPB/SARP_CSGA_merge:/input_data \
    --bind ${RKJCOLLAB}/Collabs/ortega/data/genetics/tm_r3_imp:/output_data \
    /Users/slacksa/repos/imputation_snakemake/envs/topmed_imputation.sif \
    snakemake --snakefile /repo/Snakefile \
        --configfile /proj_repo/config_topmed.yml \
        --cores 8 --until "$step"

# 1000G phase 3 v5
echo "1000G"
apptainer exec \
    --writable-tmpfs \
    --bind /Users/slacksa/repos/imputation_snakemake:/repo \
    --bind /Users/slacksa/repos/local_ancestry_analysis/imputation:/proj_repo \
    --bind ${RKJCOLLAB}/Collabs/ortega/background/from_dayam_server/Victor_Ortega_analyses_MPB/SARP_CSGA_merge:/input_data \
    --bind ${RKJCOLLAB}/Collabs/ortega/data/genetics/1000g_imp:/output_data \
    /Users/slacksa/repos/imputation_snakemake/envs/topmed_imputation.sif \
    snakemake --snakefile /repo/Snakefile \
        --configfile /proj_repo/config_1000g.yml \
        --cores 8 --until "$step"

# # HRC r1.1
echo "HRC"
apptainer exec \
    --writable-tmpfs \
    --bind /Users/slacksa/repos/imputation_snakemake:/repo \
    --bind /Users/slacksa/repos/local_ancestry_analysis/imputation:/proj_repo \
    --bind ${RKJCOLLAB}/Collabs/ortega/background/from_dayam_server/Victor_Ortega_analyses_MPB/SARP_CSGA_merge:/input_data \
    --bind ${RKJCOLLAB}/Collabs/ortega/data/genetics/hrc_imp:/output_data \
    /Users/slacksa/repos/imputation_snakemake/envs/topmed_imputation.sif \
    snakemake --snakefile /repo/Snakefile \
        --configfile /proj_repo/config_hrc.yml \
        --cores 8 --until "$step"
