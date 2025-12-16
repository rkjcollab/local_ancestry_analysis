# **Overview**

This repo contains general code for running RFMix-based local ancestry analysis
(admixture mapping and LAAA fine mapping). The repo rkjcollab/pgx_laba contains
example scripts for running the code here. Originally based on Michelle Daya's
repository mdaya/local_anc_afr_eur_on_sevenbridges.

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

# **Analysis Options**

## **1. local_anc_afr_eur (RFMix)**

This sub-folder contains all the scripts needed to calculate two-way (AFR and
EUR) local ancestry estimates using RFMix.

Bash batch scripts can be used to run each step described below. Because they
are project-specific, they are not included here. Please see examples in the
rkjcollab repo pgx_laba.

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

TO NOTE: If the process for unphased data is run on denser, imputed or WGS data
with SHAPEIT for phasing and the references included in this repo, the local
ancestry calls will likely be incorrect (for example, showing 100% of one ancestry
for all individuals). We are unsure why this happens, but know that if possible,
the pipeline should be run on WGS or imputed data as originally phased.

## **2. gwide_anc**

This sub-folder contains all the scripts needed for quality control and summary
of local ancestry estimates from RFMix after running local_anc_afr_eur.

*For admixture mapping or QC of RFMix estimates, this full section needs to*
*be run. If only want to run LAAA, then can just run first script.*

1. 1.1_calc_gwide_ancestry.sh

The first step calculates a global ancestry proportion based on the RFMix
estimates which is used as a covariate in LAAA. For this first step to be run,
have to already have results for all chromosomes from RFMix from the
local_anc_afr_eur section.

2: 1.2_calc_local_afr_ancestry.sh

*TO NOTE: Only run this step if you are planning to run admixture mapping.*

(Optional) 3: local_ancestry_qc.Rmd

This report can be run using the output from steps 1 and 2 to verify that RFMix
calls seems accurate.


## **(Optional) 3. admixture_mapping**

This sub-folder contains all the scripts needed to run admixture mapping after
the local_anc_afr_eur and gwide_anc steps have been run.

1. Make project-specific phenotype and covariate files.

*Example code for this step included in the rkjcollab/pgx_laba repo.*

2. 1.1_create_plink_input_local.sh

This step uses a study PLINK file and the RFMix results (run in batches using
the -co option) to make PLINK files that store local ancestry estimates.

3. 1.2_run_assoc_linear.sh

This step runs the admixture mapping association for each different phenotype
file provided in a given directory. Phenotype files should have the format
"{STUDY_PREFIX}_admix_map_{PHENO}_pheno(_male|_female)?.txt".

4. 1.4_calc_m_eff.sh

This step uses the PLINK files with local ancestry estimates to calculate the
numer of effective tests made in admixture mapping to inform the threshold
for statistical significance.

5. Make project-specific report summarizes results and generating peak region
files.

TO NOTE: peak region files are required for LAAA fine mapping. They must have
the file name format "{PHENO}_admixture_peak_contig_region_chr{#}(_male|_female)?.txt".

*Example code for this step included in the rkjcollab/pgx_laba repo.*


## **(Optional) 3. laaa_fine_mapping**

This sub-folder contains all the scripts needed to run local ancestry adjusted
allelic association (LAAA) fine mapping of admixture mapping peaks. Before this
step is run, the local_anc_afr_eur subfolder and step one of the gwide_anc
subfolder must be run on the input data, which should be phased data (WGS or
imputed) for fine mapping. Additionally, admixture mapping peak regions need
to be identified, either by running the admixture_mapping subfolder or by using
regions identified in another dataset.

### **To Process Phased Data (WGS or Imputed)**

1. Make project-specific phenotype file.

The file should be tab-delimited, have a .txt extension, and include column names.

*Example code for this step included in the rkjcollab/pgx_laba repo.*

2. 1.1_create_dose_frames.sh

This script uses the admixture mapping peak region file
“{PHENO}_admixture_peak_contig_region_chr{#}.txt” (described above) to
make dose frame files summarizing the RFMix local ancestry in peak regions.

3. 1.2_run_models.sh

The LAAA model is run using allele, ancestry, and allele-ancestry dose
frames are made.

4. Make project-specific report summarizing results.

*Example code for this step included in the rkjcollab/pgx_laba repo.*

# **Additional Data**

### shapeit_input.tar.gz

Present in repository and automatically loaded into container. This is a
HapMap-based genetic map used for SHAPEIT phasing. Originally from Michelle
Daya's GitHub repository.

### rfmix_input.tar.gz

This file is too large to keep on repository (~2 GB). Please contact sdslack
for access if needed.

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

In addition to the above, the rfmix_input.tar.gz also includes file
igsr-1000_genomes_phase_3_release.tsv, which provides phenotype information
about the reference population used.
