# **Setup**

Create the Apptainer/Singularity container using "local_ancestry_analysis.def":

``` bash
apptainer build local_ancestry_analysis.sif local_ancestry_analysis.def

```

Note that the container was built for a Linux AMD64 (x86_64) system.

# **Steps**

## **local_anc_afr_eur**

This sub-folder contains all the scripts needed to calculated two-way (AFR and
EUR) local ancestry estimates using RFMix.

Bash scripts and steps are described below. For each bash script, there is a
test batch script with an example of how it can be run. These are located under
batch_files/job_<name_of_script>.batch. They are largely in SLURM batch
submission format, but can be run as simple bash scripts.

### **To Process Phased Data (WGS or Imputed)**

1. prep_input_phased.sh

2. run_hg38_phased_conversion.sh

3. run_rfmix.sh

Before running the first time, need to unzip rfmix_input.tar.gz (described in
**Additional Data** secion):

```bash
tar -zxvf rfmix_input.tar.gz 

```


### **To Process Unphased Data (Chip)**

1. prep_input_unphased.sh

2. run_shapeit.sh

3. run_rfmix.sh


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
