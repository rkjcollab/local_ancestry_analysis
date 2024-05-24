#!/bin/bash
#SBATCH --partition bigmem
#SBATCH --time=1440
#SBATCH --ntasks=1 
#SBATCH --nodes=1 
#SBATCH --mem=1024000 
#SBATCH --mail-user=michelle.daya@ucdenver.edu
#SBATCH --account=barnescaapa-dayam

#Set parameters
if [ "$#" -eq  "0" ]
then
    echo "Usage: ${0##*/} <chr>"
    exit
fi
chr=$1


#Make sure output directories exist
mkdir ../data/output
mkdir ../data/output/rfmix
rm -r ../data/output/rfmix/chr${chr}
mkdir ../data/output/rfmix/chr${chr}
#Move SNPs file
cp tmp_chr${chr}_snps_keep.txt ../data/output/rfmix/chr${chr}/snps.txt

c=1
b=1
while [  $c -lt 961 ]; do
   let begin_haplo=c
   let end_haplo=c+19
   bash run_rfmix_popphased_rosie.sh $chr $begin_haplo $end_haplo $b
   let c=c+20
   let b=b+1
done
let begin_haplo=c
let end_haplo=982
bash run_rfmix_popphased_rosie.sh $chr $begin_haplo $end_haplo $b

cd ../data/output/rfmix/chr${chr}
paste local_ancestry_popphased_batch1.0.Viterbi.txt \
   local_ancestry_popphased_batch2.0.Viterbi.txt \
   local_ancestry_popphased_batch3.0.Viterbi.txt \
   local_ancestry_popphased_batch4.0.Viterbi.txt \
   local_ancestry_popphased_batch5.0.Viterbi.txt \
   local_ancestry_popphased_batch6.0.Viterbi.txt \
   local_ancestry_popphased_batch7.0.Viterbi.txt \
   local_ancestry_popphased_batch8.0.Viterbi.txt \
   local_ancestry_popphased_batch9.0.Viterbi.txt \
   local_ancestry_popphased_batch10.0.Viterbi.txt \
   local_ancestry_popphased_batch11.0.Viterbi.txt \
   local_ancestry_popphased_batch12.0.Viterbi.txt \
   local_ancestry_popphased_batch13.0.Viterbi.txt \
   local_ancestry_popphased_batch14.0.Viterbi.txt \
   local_ancestry_popphased_batch15.0.Viterbi.txt \
   local_ancestry_popphased_batch16.0.Viterbi.txt \
   local_ancestry_popphased_batch17.0.Viterbi.txt \
   local_ancestry_popphased_batch18.0.Viterbi.txt \
   local_ancestry_popphased_batch19.0.Viterbi.txt \
   local_ancestry_popphased_batch20.0.Viterbi.txt \
   local_ancestry_popphased_batch21.0.Viterbi.txt \
   local_ancestry_popphased_batch22.0.Viterbi.txt \
   local_ancestry_popphased_batch23.0.Viterbi.txt \
   local_ancestry_popphased_batch24.0.Viterbi.txt \
   local_ancestry_popphased_batch25.0.Viterbi.txt \
   local_ancestry_popphased_batch26.0.Viterbi.txt \
   local_ancestry_popphased_batch27.0.Viterbi.txt \
   local_ancestry_popphased_batch28.0.Viterbi.txt \
   local_ancestry_popphased_batch29.0.Viterbi.txt \
   local_ancestry_popphased_batch30.0.Viterbi.txt \
   local_ancestry_popphased_batch31.0.Viterbi.txt \
   local_ancestry_popphased_batch32.0.Viterbi.txt \
   local_ancestry_popphased_batch33.0.Viterbi.txt \
   local_ancestry_popphased_batch34.0.Viterbi.txt \
   local_ancestry_popphased_batch35.0.Viterbi.txt \
   local_ancestry_popphased_batch36.0.Viterbi.txt \
   local_ancestry_popphased_batch37.0.Viterbi.txt \
   local_ancestry_popphased_batch38.0.Viterbi.txt \
   local_ancestry_popphased_batch39.0.Viterbi.txt \
   local_ancestry_popphased_batch40.0.Viterbi.txt \
   local_ancestry_popphased_batch41.0.Viterbi.txt \
   local_ancestry_popphased_batch42.0.Viterbi.txt \
   local_ancestry_popphased_batch43.0.Viterbi.txt \
   local_ancestry_popphased_batch44.0.Viterbi.txt \
   local_ancestry_popphased_batch45.0.Viterbi.txt \
   local_ancestry_popphased_batch46.0.Viterbi.txt \
   local_ancestry_popphased_batch47.0.Viterbi.txt \
   local_ancestry_popphased_batch48.0.Viterbi.txt \
   local_ancestry_popphased_batch49.0.Viterbi.txt \
   -d' ' | sed 's/  / /g' > local_ancestry.0.Viterbi.txt

paste local_ancestry_popphased_batch1.allelesRephased0.txt \
   local_ancestry_popphased_batch2.allelesRephased0.txt \
   local_ancestry_popphased_batch3.allelesRephased0.txt \
   local_ancestry_popphased_batch4.allelesRephased0.txt \
   local_ancestry_popphased_batch5.allelesRephased0.txt \
   local_ancestry_popphased_batch6.allelesRephased0.txt \
   local_ancestry_popphased_batch7.allelesRephased0.txt \
   local_ancestry_popphased_batch8.allelesRephased0.txt \
   local_ancestry_popphased_batch9.allelesRephased0.txt \
   local_ancestry_popphased_batch10.allelesRephased0.txt \
   local_ancestry_popphased_batch11.allelesRephased0.txt \
   local_ancestry_popphased_batch12.allelesRephased0.txt \
   local_ancestry_popphased_batch13.allelesRephased0.txt \
   local_ancestry_popphased_batch14.allelesRephased0.txt \
   local_ancestry_popphased_batch15.allelesRephased0.txt \
   local_ancestry_popphased_batch16.allelesRephased0.txt \
   local_ancestry_popphased_batch17.allelesRephased0.txt \
   local_ancestry_popphased_batch18.allelesRephased0.txt \
   local_ancestry_popphased_batch19.allelesRephased0.txt \
   local_ancestry_popphased_batch20.allelesRephased0.txt \
   local_ancestry_popphased_batch21.allelesRephased0.txt \
   local_ancestry_popphased_batch22.allelesRephased0.txt \
   local_ancestry_popphased_batch23.allelesRephased0.txt \
   local_ancestry_popphased_batch24.allelesRephased0.txt \
   local_ancestry_popphased_batch25.allelesRephased0.txt \
   local_ancestry_popphased_batch26.allelesRephased0.txt \
   local_ancestry_popphased_batch27.allelesRephased0.txt \
   local_ancestry_popphased_batch28.allelesRephased0.txt \
   local_ancestry_popphased_batch29.allelesRephased0.txt \
   local_ancestry_popphased_batch30.allelesRephased0.txt \
   local_ancestry_popphased_batch31.allelesRephased0.txt \
   local_ancestry_popphased_batch32.allelesRephased0.txt \
   local_ancestry_popphased_batch33.allelesRephased0.txt \
   local_ancestry_popphased_batch34.allelesRephased0.txt \
   local_ancestry_popphased_batch35.allelesRephased0.txt \
   local_ancestry_popphased_batch36.allelesRephased0.txt \
   local_ancestry_popphased_batch37.allelesRephased0.txt \
   local_ancestry_popphased_batch38.allelesRephased0.txt \
   local_ancestry_popphased_batch39.allelesRephased0.txt \
   local_ancestry_popphased_batch40.allelesRephased0.txt \
   local_ancestry_popphased_batch41.allelesRephased0.txt \
   local_ancestry_popphased_batch42.allelesRephased0.txt \
   local_ancestry_popphased_batch43.allelesRephased0.txt \
   local_ancestry_popphased_batch44.allelesRephased0.txt \
   local_ancestry_popphased_batch45.allelesRephased0.txt \
   local_ancestry_popphased_batch46.allelesRephased0.txt \
   local_ancestry_popphased_batch47.allelesRephased0.txt \
   local_ancestry_popphased_batch48.allelesRephased0.txt \
   local_ancestry_popphased_batch49.allelesRephased0.txt \
   -d' ' | sed 's/ //g' > local_ancestry.allelesRephased0.txt

