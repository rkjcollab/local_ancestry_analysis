---
title: "README"
---

Pipeline for local ancestry inference of an African and European admixed population, for which the TGP YRI and CEU populations are good ancestral population proxies. Note that this pipeline assumes an hg19 build of the genome. 

## Software and other requirements

- ShapeIT (v2.r837 was last used)
- RFMix (v1.5.4 was last used)
- VCFTools (v0.1.13 was used)
- ADMIXTURE (Version 1.3.0 was last used)
- R Studio setup to build R markdown files
- TGP and map input files should be present in the CAAPA local ancestry inference (/gpfs/share/barnescaapa/caapa_local_ancestry/data/input)

### Details of TGP and map input files

To save uncessary storage space and copying/soft link files, several input files created for the CAAPA local ancestry inference are used from the CAAPA local ancestry inference folder. They were initially created as part of that project, as follows:

### genetic_map_hapmap

chr\<nr\>.txt files, required for ShapeIt, as per the following ShapeIt page:
https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html#gmap .

The hapmap genetic map file was downloaded from: http://www.shapeit.fr/files/genetic_map_b37.tar.gz on 29 April 2016.

### genetic_map_tgp

Genetic map files that are required by RFMix. 

The map files used were downloaded from https://github.com/joepickrell/1000-genomes-genetic-maps/tree/master/interpolated_from_hapmap (downloaded oned 29 April 2016).
The files were renamed to chr1.txt, chr2.txt, ..., chr22.txt

### tgp

Reference (ancestral) populations required by RFmix - 99 CEU and 108 YRI. 

Downloaded the VCF files in ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ to data/raw/tgp_release_20130502. Created the following file from the first worksheet of http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/20130606_sample_info.xlsx:
<i>ceu_yri_ids.txt</i> in this same directory, containing a list of all CEU and YRI IDs. Ran the script <code>create_tgp_input_files.sh</code> in the scripts directory, to extract the already phased TGP CEU and YRI subjects into IMPUTE file format (which is easy to later on merge with the admixed ShapeIT phased files).

## Step 1: Setup

Run <code>1_setup.sh</code> on Rosie. 

Copy (or create a soft link) to the QC version of the GWAS array PLINK files to data/input/admixed.fam (.bim .bed).

Note that ShapeIT does not allow for SNPs with duplicate positions, so remove
these first.

## Step 2: Create admixed phased files by chromosome

(Note: ths is a slurm script)

<code>2_run_shapeit.sh \<chr\></code>

## Step 3: Run RFMix by chromosome 

This slrum script also takes care of creating the RFMix input files.

<code>3_run_rfmix.sh \<chr\></code>

## Step 4: Estimate genome-wide ancestry 

Estimate genome-wide ancestry from RFMix output, and estimate ADMIXTURE genome-wide ancestry.
Merge the genome-wide ancestries.
This is all done by running the following slurm script:

<code>4_calc_gwide_ancestry.sh</code>

## Step 5: Calculate the proportion of African ancestry at each local ancestry segment

(Note: this is NOT a slurm script)

<code>5_calc_local_afr_ancestry.sh</code>

## Step 6: Do local ancestry QC

Copy over the following files from Rosie to your local machine, maintaining the directory structure:

* data/output/merged_gwide.txt
* data/output/local_ancestry_segments.txt
* data/output/rfmix/chr*/snps.txt (easiest is to zip the entire rfmix directory - <code>zip -r rfmix.zip rfmix</code> - and copy that over and then extract it)

Run <i>6_local_ancestry_qc.Rmd</i> in R Studio


