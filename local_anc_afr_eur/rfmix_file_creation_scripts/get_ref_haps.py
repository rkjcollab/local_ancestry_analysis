import os, sys

snp_keep_file = open(sys.argv[1])
hap_file = out_file = open(sys.argv[2])
out_file = open(sys.argv[3], "w")
hap_file.readline() #skip the header line

prev_snp = ""
hap_pos = ""
for snp_line in snp_keep_file:
    e=snp_line.strip().split()
    snp_pos = e[0]
    while hap_pos != snp_pos:
        hap_line=hap_file.readline()
        e = hap_line.strip().split()
        snp = e[0]
        hap_pos = snp.split(":")[1]
        if (snp != prev_snp) and (hap_pos == snp_pos):
            for i in range(1, len(e)-1):
                out_file.write(e[i] + " ")
            out_file.write(e[len(e)-1] + "\n")
        prev_snp = snp

snp_keep_file.close()
hap_file.close()
out_file.close()
