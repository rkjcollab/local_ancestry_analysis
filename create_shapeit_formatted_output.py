import sys
import gzip

chr=str(sys.argv[1])
imputed_file_name="../../csga_imputation/data/output/chr" + chr + ".dose.vcf.gz"
sample_out_file_name="../data/intermediate/shapeit_output/chr" + chr + ".samples"
hap_out_file_name="../data/intermediate/shapeit_output/chr" + chr + ".haps"
bed_file_name="tmp_chr" + chr + "_out.bed"

#Process the imputed file
sample_out_file=open(sample_out_file_name, "w")
f=gzip.open(imputed_file_name, 'rb')
chr_str="chr" + chr
vcf_dose_map={}
for line in f:
    el=line.strip().split()
    #Write sample IDs to sample file
    if el[0] == "#CHROM": 
        for i in range(9,len(el)):
            sample_out_file.write(el[i] + "\n")
        sample_out_file.close()
    if el[0] == chr_str:
        pos = el[1]
        variant_id = el[2]
        ref = el[3]
        alt = el[4]
        allele_str = ""
        for i in range(9,len(el)):
            geno=el[i].strip().split(":")[0]
            alleles=geno.split("|")
            allele_str = allele_str + " " + alleles[0] + " " + alleles[1]
        vcf_dose_map[pos] = [variant_id, ref, alt, allele_str]
f.close()

#Process the cross map output file
f=open(bed_file_name)
hap_out_file=open(hap_out_file_name, "w")
for line in f:
    el=line.strip().split()
    hg38_pos=el[3].split(":")[1]
    hg19_pos=el[1]
    if hg38_pos in vcf_dose_map.keys():
        val = vcf_dose_map[hg38_pos]
        variant_id = val[0]
        pos = hg19_pos
        ref = val[1]
        alt = val[2]
        allele_str = val[3]
        hap_out_file.write(chr + " " + variant_id + " " + pos + " " + ref + " " + alt + allele_str + "\n")
f.close()
hap_out_file.close()
