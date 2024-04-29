#!/bin/bash

#Set parameters
if [ "$#" -eq  "0" ]
then
    echo "Usage: ${0##*/} <chr>"
    exit
fi
chr=$1
begin_haplo=$2
end_haplo=$3
batch=$4

cd ../data/input/rfmix/chr${chr}
cut -f${begin_haplo}-${end_haplo},983-1396 -d' ' classes.txt > classes_batch${batch}.txt
cut -c${begin_haplo}-${end_haplo},983-1396 alleles.txt > alleles_batch${batch}.txt
cd ../../../../scripts

#Run RFMix
RFMix_PopPhased -a ../data/input/rfmix/chr${chr}/alleles_batch${batch}.txt \
                  -p ../data/input/rfmix/chr${chr}/classes_batch${batch}.txt \
                  -m ../data/input/rfmix/chr${chr}/snp_locations.txt \
                  -co 1 -o ../data/output/rfmix/chr${chr}/local_ancestry_popphased_batch${batch}

