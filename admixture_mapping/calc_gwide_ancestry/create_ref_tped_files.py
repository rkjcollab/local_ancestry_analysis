import sys

chr=sys.argv[1]
admix_input_dir=sys.argv[2]
rfmix_working_dir=sys.argv[3]

snp_file=open(admix_input_dir + "/chr" + chr + ".keep")
snp_hash={}
for line in snp_file:
    snp=line.strip()
    snp_hash[snp] = ""
snp_file.close()

hap_file=open(rfmix_working_dir + "/chr" + chr + ".hap")
out_file=open(rfmix_working_dir + "/chr" + chr + ".tped", "w")
for line in hap_file:
    e=line.strip().split()
    snp=e[0]
    pos=e[1]
    a1=e[2]
    a2=e[3]
    if pos in snp_hash:
        out_file.write(chr + " " + snp + " 0 " + pos)
        for i in range(4, len(e)):
            a=e[i]
            if a == "0":
                out_file.write(" " + a1)
            else:
                out_file.write(" " + a2)
        out_file.write("\n")
hap_file.close()
out_file.close()
