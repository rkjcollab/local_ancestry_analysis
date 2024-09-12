#!/bin/bash

if [ "$#" -eq  "0" ]
then
   echo "Usage: ${0##*/} <path_to_vcf_input> <output_dir_prefix> <nr_threads>"
   echo "      <nr_threads> <code_dir> <min_rsq> <genome_build_of_input>"
   echo "Script converts WGS or imputed input to SHAPEIT-formatted output,"
   echo "maintaining phasing for local ancestry pipeline. Input"
   echo "should be a single-chromosome VCF file in either hg19 or hg38."
   echo "<genome_build_of_input> must either be 'hg19' or 'hg38'. If input"
   echo "is WGS data, set <min_rsq> to 0."
   exit
fi

vcf_input=$1
out_dir=$2
nr_threads=$3
code_dir=$4
min_rsq=$5  # if processing WGS, set to 0
build=$6  # either hg38 or hg19

# Make sure genome build only hg19 or hg38
if [[ "$build" != "hg19" && "$build" != "hg38" ]]; then
	echo "Input build must be either 'hg38' or 'hg19'."
	exit
fi

file_id=`basename $vcf_input | sed 's/.vcf.gz//'`

# Apply the Rsq filter (if zero all the SNPs will be retained), so 0
# should be used for WGS input.
if [ "$min_rsq" == "0" ] 
then
   echo "<min_rsq> set to 0, likely for WGS data"
   touch ${out_dir}/tmp_${file_id}_r2_failed_variants.txt
   # create an empty file 
else
   echo -e 'CHROM\tPOS\tREF\tALT\tR2' > ${out_dir}/tmp_${file_id}_r2.txt
      # needs to be single ''
   bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%R2\n' $vcf_input >> \
      ${out_dir}/tmp_${file_id}_r2.txt

   cat ${code_dir}/shapeit_formatting_scripts/filter_rsq.R | R --vanilla --args \
      ${out_dir}/tmp_${file_id}_r2.txt \
      ${out_dir}/tmp_${file_id}_r2_failed_variants.txt \
      $min_rsq
fi
if [ -s "${out_dir}/tmp_${file_id}_r2_failed_variants.txt" ]
then
   rsq_vcf_file=${out_dir}/tmp_${file_id}_min_maf_with_snp_id_min_rsq.vcf.gz
   bcftools view -e "ID=@${out_dir}/tmp_${file_id}_r2_failed_variants.txt" \
      $vcf_input -Oz -o $rsq_vcf_file
      # needs to be double ""
else
   echo "<min_rsq> set to 0, likely for WGS data"
   rsq_vcf_file=$vcf_input
fi

# If input is hg38 liftover to hg19, based on code from topmed_imputation repo
if [ "$build" == "hg38" ]; then
    echo "Input build is hg38. Lifting over to hg19."
    # Get lifted over coordinates of variants 
    zcat $rsq_vcf_file | grep -v ^# | cut -f1 > "${out_dir}/tmp_${file_id}_chr_col.txt"
    zcat $rsq_vcf_file | grep -v ^# | cut -f2 > "${out_dir}/tmp_${file_id}_pos_col.txt"
    zcat $rsq_vcf_file | grep -v ^# | cut -f3 > "${out_dir}/tmp_${file_id}_snp_col.txt"
    chr=`head -1 ${out_dir}/tmp_${file_id}_chr_col.txt | sed 's/chr//'`
    paste "${out_dir}/tmp_${file_id}_chr_col.txt" \
      "${out_dir}/tmp_${file_id}_pos_col.txt" "${out_dir}/tmp_${file_id}_pos_col.txt" \
      "${out_dir}/tmp_${file_id}_snp_col.txt" > "${out_dir}/tmp_chr${chr}_in.bed"

    CrossMap bed "${code_dir}/shapeit_formatting_scripts/hg38ToHg19.over.chain" \
                "${out_dir}/tmp_chr${chr}_in.bed"  \
                "${out_dir}/tmp_chr${chr}_out.bed"

    # Create bed file to annotate the hg19 position
    cat "${code_dir}/shapeit_formatting_scripts/create_hg19_annot_bed.R" | \
        R --vanilla --args $chr \
        "${out_dir}/tmp_chr${chr}_out.bed" \
        "${out_dir}/tmp_chr${chr}_hg19_annot.txt"
    bgzip "${out_dir}/tmp_chr${chr}_hg19_annot.txt"
    tabix -s1 -b2 -e3 "${out_dir}/tmp_chr${chr}_hg19_annot.txt.gz"

    # Annotate the hg19 column
    hg19_annot_vcf_file="${out_dir}/tmp_${file_id}_min_maf_with_snp_id_min_rsq_hg19_annot.vcf"
    bcftools annotate \
      -a "${out_dir}/tmp_chr${chr}_hg19_annot.txt.gz" \
      -c CHROM,FROM,TO,REF,ALT,HG19 \
      -h <(echo '##INFO=<ID=HG19,Number=1,Type=Integer,Description="hg19 position">') \
      $rsq_vcf_file > $hg19_annot_vcf_file

    # Filter out variants without an hg19 annotation (line ends with .)
    bcftools query -f '%ID\t%INFO/HG19\n' $hg19_annot_vcf_file > \
        "${out_dir}/tmp_${file_id}_hg19.txt"
    grep -v "\.$" "${out_dir}/tmp_${file_id}_hg19.txt" |
        cut -f1 > \
        "${out_dir}/tmp_${file_id}_hg19_variants.txt"
    hg19_filtered_vcf_file="${out_dir}/tmp_${file_id}_hg19_only.vcf"
    bcftools view -i ID=@"${out_dir}/tmp_${file_id}_hg19_variants.txt" \
        $hg19_annot_vcf_file \
        -Ov -o $hg19_filtered_vcf_file

    # Update POS column with hg19 position
    hg19_pos_vcf_file="${out_dir}/tmp_${file_id}_hg19_pos.vcf"
    grep "^#" $hg19_filtered_vcf_file > $hg19_pos_vcf_file
    grep -v "^#" ${hg19_filtered_vcf_file} | cut -f1 > ${hg19_filtered_vcf_file}.c1 #chr
    bcftools query -f "%INFO/HG19\n" $hg19_filtered_vcf_file > \
        ${hg19_filtered_vcf_file}.c2 #pos

    grep -v "^#" $hg19_filtered_vcf_file | cut -f 3- >  ${hg19_filtered_vcf_file}.c3_onwards
    paste ${hg19_filtered_vcf_file}.c1 \
        ${hg19_filtered_vcf_file}.c2 \
        ${hg19_filtered_vcf_file}.c3_onwards >> $hg19_pos_vcf_file

else
	echo "Input build is already hg19."

    hg19_pos_vcf_file="${out_dir}/tmp_${file_id}_hg19_pos.vcf"
    cp $rsq_vcf_file $hg19_pos_vcf_file
fi

# Sort by position
hg19_sorted_vcf_file="${out_dir}/tmp_${file_id}_hg19_sorted.vcf"
bcftools sort $hg19_pos_vcf_file -o $hg19_sorted_vcf_file

# Convert the hg19 VCF to a shapeit format haps file chr${chr}.haps
grep ^# -v $hg19_sorted_vcf_file > ${hg19_sorted_vcf_file}.genos
cut -f1 ${hg19_sorted_vcf_file}.genos > ${hg19_sorted_vcf_file}.haps.c1 #chr
cut -f3 ${hg19_sorted_vcf_file}.genos > ${hg19_sorted_vcf_file}.haps.c2 #snp_id
cut -f2 ${hg19_sorted_vcf_file}.genos > ${hg19_sorted_vcf_file}.haps.c3 #pos
cut -f4 ${hg19_sorted_vcf_file}.genos > ${hg19_sorted_vcf_file}.haps.c4 #ref
cut -f5 ${hg19_sorted_vcf_file}.genos > ${hg19_sorted_vcf_file}.haps.c5 #alt
cut -f 10- ${hg19_sorted_vcf_file}.genos > ${hg19_sorted_vcf_file}.haps.c6_onwards #genotypes

paste -d' ' ${hg19_sorted_vcf_file}.haps.c1 \
    ${hg19_sorted_vcf_file}.haps.c2 \
    ${hg19_sorted_vcf_file}.haps.c3 \
    ${hg19_sorted_vcf_file}.haps.c4 \
    ${hg19_sorted_vcf_file}.haps.c5 \
    ${hg19_sorted_vcf_file}.haps.c6_onwards > "${out_dir}/chr${chr}.haps"

# Create the samples file chr${chr}.samples
sample_line=`grep \#CHROM $hg19_sorted_vcf_file | cut -f 10-` 
echo "ID2" > "${out_dir}/chr${chr}.samples"
echo "0" >> "${out_dir}/chr${chr}.samples"
for value in $sample_line
do
    echo $value >> "${out_dir}/chr${chr}.samples"
done

# Clean up
rm ${out_dir}/tmp_*
