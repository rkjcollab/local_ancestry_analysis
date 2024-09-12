# **Setup**

*Note that these steps should only need to be run once.*

Create the Apptainer/Singularity container using "local_ancestry_analysis.def":

``` bash
apptainer build local_ancestry_analysis.sif local_ancestry_analysis.def

```

Note that the container was built for a Linux AMD64 (x86_64) system.

Unzip rfmix_input.tar.gz (described in **Additional Data** secion):

```bash
tar -zxvf rfmix_input.tar.gz

```

If move this input to same local folder as repository, path to its location will
not need to be updated in scripts. This is added to the .gitignore, so it can't
be committed or pushed to the repository.

```bash
mv -r rfmix_input/ /path/to/repo

```

# **Steps**

## **local_anc_afr_eur**

This sub-folder contains all the scripts needed to calculated two-way (AFR and
EUR) local ancestry estimates using RFMix.

Bash scripts and steps are described below. For each bash script, there is a
test batch script with an example of how it can be run. These are located under
batch_files/job_<name_of_script>.batch. They are largely in SLURM batch
submission format, but can be run as simple bash scripts.

Many additional details are included in the bash and batch scripts.

### **To Process Phased Data (WGS or Imputed)**

This is likely the desired input data type for ultimately running LAAA fine
mapping.

1. prep_input_phased.sh

2. run_hg38_phased_conversion.sh

3. run_rfmix.sh

RFMix collapse option should be set to -o for uncollapsed output (used for LAAA
input) and set to -co for collapsed output (used for Admixture Mapping inout).

### **To Process Unphased Data (Chip)**

This is likely the desired input data type for ultimately running admixture
mapping.

1. prep_input_unphased.sh

2. run_shapeit.sh

3. run_rfmix.sh

RFMix collapse option should be set to -o for uncollapsed output (used for LAAA
fine mapping input) and set to -co for collapsed output (used for admixture
mapping inout).

## **admixture_mapping**

*TODO: to be added*


## **laaa_fine_mapping**

*TODO: to be added*


# **Additional Data**

### shapeit_input.tar.gx

Present in repository and automatically loaded into container.

*TODO: add info here!*

### rfmix_input.tar.gz

Too large to keep on repository.

*TODO: find better way tot share?*
*TODO: add details!*
