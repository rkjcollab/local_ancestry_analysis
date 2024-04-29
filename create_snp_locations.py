import os, sys

chr=sys.argv[1]

snp_keep_file = open("tmp_chr" + chr + "_snps_keep.txt")
snps_hash = {}
for line in snp_keep_file:
    e=line.strip().split()
    snps_hash[e[0]] = ""
snp_keep_file.close()

map_file = open("/gpfs/share/barnescaapa/caapa_local_ancestry/data/input/genetic_map_tgp/chr" + chr + ".txt")
cm_hash = {}
for line in map_file:
    e=line.strip().split()
    pos=e[1]
    cm=e[2]
    if pos in snps_hash.keys():
        cm_hash[pos] = cm
map_file.close()

hap_file = open("../data/intermediate/shapeit_output/chr" + chr + ".haps")
out_file = open("../data/input/rfmix/chr" + chr + "/snp_locations.txt", "w")
i = 0
for line in hap_file:
    e=line.strip().split()
    pos = e[2]
    if pos in snps_hash.keys():
        if pos in cm_hash:
            out_file.write(cm_hash[pos] + "\n")
    i = i + 1
hap_file.close()
out_file.close()
