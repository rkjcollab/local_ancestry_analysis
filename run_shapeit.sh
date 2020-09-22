#!/bin/bash

bed_fn=$1
genetic_map_input_dir=$2
nr_threads=$3

#Set bim and fam file name
bim_fn=`echo $bed_fn | sed 's/.bed/.bim/'`
fam_fn=`echo $bed_fn | sed 's/.bed/.fam/'`

#Set chr variable 
#use the first column of the first line of the bim file to get this
head -1 $bim_fn | awk '{print $1}' | sed 's/chr//'

#Run shapeit
shapeit --input-bed $bed_fn $bim_fn $fam_fn \
        --input-map ${genetic_map_input_dir}/chr${chr}.txt \
        --thread $nr_threads \
        --output-max chr${chr}.haps chr${chr}.samples
