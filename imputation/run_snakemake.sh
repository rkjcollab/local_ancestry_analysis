#!/bin/bash

#SBATCH --nodes=1
#SBATCH --partition=amilan
#SBATCH --ntasks=22
#SBATCH --job-name="hrc_imp"  # used as part of output name
#SBATCH -o "/scratch/alpine/sslack@xsede.org/ortega/imputation/%x/imputed/job_snakemake.out"
#SBATCH -e "/scratch/alpine/sslack@xsede.org/ortega/imputation/%x/imputed/job_snakemake.err"
#SBATCH --account=amc-general
#SBATCH --time=02:00:00
#SBATCH --mem=60G  # should be size of largest N zip files processing at same time
#SBATCH --qos=normal

# This is currently called inside of the limactl amd64.
# Set paths and options in config.yml file. Paths are automatically 
# relative to snakefile location.

# TO NOTE: can add --dry-run to make sure not re-running additional steps

# Uncomment step want to run through
# step="submit_initial_input"  # pre-imp QC
# step="submit_fix_strands"  # submit imp
step="concat_convert_to_plink"  # unzip, clean, & merge

# Set number cores or get from SLURM
n_cores="$SLURM_NTASKS"

# Set base dirs
# Uncomment if local
# RKJCOLLAB="/Users/slacksa/Library/CloudStorage/OneDrive-TheUniversityofColoradoDenver"
# base_data="${RKJCOLLAB}/Collabs/ortega"
# base_code="/Users/slacksa/repos"

# Uncomment if Alpine
base_data="/scratch/alpine/sslack@xsede.org/ortega/imputation"
base_code="/projects/sslack@xsede.org/repos"

# Run snakemake
# --bind ${base_data}/background/from_dayam_server/Victor_Ortega_analyses_MPB/SARP_CSGA_merge:/input_data \
apptainer exec \
    --writable-tmpfs \
    --bind ${base_code}/imputation_snakemake:/repo \
    --bind ${base_code}/local_ancestry_analysis/imputation:/proj_repo \
    --bind ${base_data}/${SLURM_JOB_NAME}:/output_data \
    ${base_code}/imputation_snakemake/envs/topmed_imputation.sif \
    snakemake --snakefile /repo/Snakefile \
        --configfile /proj_repo/config_${SLURM_JOB_NAME}.yml \
        --cores "$n_cores" --until "$step"
