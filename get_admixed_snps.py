import os, sys

ref_fn=sys.argv[1]
snp_fn=sys.argv[2]
out_fn=sys.argv[3]

#Define a function to flip alleles
def flip (a):
    f = "0"
    if a == "A":
        f = "T"
    if a == "T":
        f = "A"
    if a == "C":
        f = "G"
    if a == "G":
        f = "C"
    return (f)

#Create SNP map file
ref_file = open(ref_fn)
snp_map={} #key:snp_pos, val=a1a2
for line in ref_file:
    e = line.strip().split()
    key = e[0]
    val = e[2] + e[3]
    snp_map[key] = val
ref_file.close()

#Read the admixed phased file and determine which SNPs to keep
snp_file = open(snp_fn)
out_file = open(out_fn, "w")
for line in snp_file:
    e = line.strip().split()
    snp = e[0]
    key = e[1]
    val = e[2] + e[3]
    val_swapped = e[3] + e[2]
    flipped_val = flip(e[2]) + flip(e[3])
    flipped_val_swapped = flip(e[3]) + flip(e[2])
    if (key in snp_map) and (val != "AT") and (val != "TA") and (val != "CG") and val != ("GC"):
        ref_val = snp_map[key]
        if (val == ref_val) or (flipped_val == ref_val):
            out_file.write(key + " 0 " + ref_val[0] + " " + ref_val[1] + "\n")
        elif (val_swapped == ref_val) or (flipped_val_swapped == ref_val):
            out_file.write(key + " 1 " + ref_val[0] + " " + ref_val[1] + "\n")
out_file.close()
