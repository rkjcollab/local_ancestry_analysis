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

This sub-folder contains all the scripts needed to calculate two-way (AFR and
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

This sub-folder contains all the scripts needed to run local ancestry adjusted
allelic association (LAAA) fine mapping of admixture mapping peaks. Before this
step is run, the local_anc_afr_eur subfolder must be run on the input data,
which should be phased data (WGS or imputed) for fine mapping. Additionally,
admixture mapping peak regions need to be identified, either by running the
admixture_mapping subfolder or by using regions identified in another dataset.

### **To Process Phased Data (WGS or Imputed)**

1. 1.1_calc_rfmix_gwide_ancestry.sh

The first step calculates a global ancestry proportion based on the RFMix
estimates which is used as a covariate in LAAA. For this first step to be run,
you will have to already have results for all chromosomes from RFMix.

*TODO: add discussion here for alternative if already have global ancestry.*

2. 2.1_create_dose_frames.sh

The second step will use the peak region files from admixture mapping (file name
format “{PHENO}_admixture_peak_contig_region_chr{#}.txt”) as inputs.

*3. 2.2_make_laaa_pheno.R*

*TO NOTE: this script is not flexible.* The third step is to make the phenotype
file to be used for analysis. The script provided in the repo is an example 
script used to create a file in the correct format. The file should be tab-
delimited, have a .txt extension, and include column names.

4. 2.3_run_models.sh

The fourth step will run the LAAA model.

# **Additional Data**

### shapeit_input.tar.gz

Present in repository and automatically loaded into container.

*TODO: add info here!*

### rfmix_input.tar.gz

Too large to keep on repository.

*TODO: find better way to share?*

The files in rfmix_input folder were prepared by Michelle Daya, and these notes
were included in her README:

<em>
The data is from the following sources:

+ genetic_map_hapmap

chr\<nr\>.txt files, required for ShapeIt, as per the following ShapeIt page:
https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html#gmap .

The hapmap genetic map file was downloaded from:
	http://www.shapeit.fr/files/genetic_map_b37.tar.gz on 29 April 2016.

+ genetic_map_tgp

Genetic map files that are required by RFMix. 

The map files used were downloaded from:
	 https://github.com/joepickrell/1000-genomes-genetic-maps/tree/master/interpolated_from_hapmap on 29 April 2016.

The files were renamed to chr1.txt, chr2.txt, ..., chr22.txt

+ tgp

Reference (ancestral) populations required by RFmix - 99 CEU and 108 YRI. 

Downloaded the VCF files in ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/
to data/raw/tgp_release_20130502. Created the following file from the first
worksheet of http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_sample_info.xlsx:
<i>ceu_yri_ids.txt</i> in this same directory, containing a list of all CEU and
YRI IDs. Ran the script <code>create_tgp_input_files.sh</code> in the scripts
directory, to extract the already phased TGP CEU and YRI subjects into IMPUTE
file format (which is easy to later on merge with the admixed ShapeIT phased
files).
</em>