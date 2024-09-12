import os, sys

snp_keep_file = open(sys.argv[1])
hap_file = open(sys.argv[2])
out_file = open(sys.argv[3], "w")

hap_pos = ""
for snp_line in snp_keep_file:
    e=snp_line.strip().split()
    snp_pos=e[0]
    flip_code=e[1]
    while hap_pos != snp_pos:
        hap_line=hap_file.readline()
        e=hap_line.strip().split()
        hap_pos = e[2]
        if hap_pos == snp_pos:
            if flip_code == "0":
                for i in range(5, len(e)-1):
                    out_file.write(e[i] + " ")
                out_file.write(e[len(e)-1] + "\n")
            else:
                for i in range(5, len(e)-1):
                    if e[i] == "0":
                        out_file.write("1 ")
                    else:
                        out_file.write("0 ")
                if e[len(e)-1] == "0":
                    out_file.write("1\n")
                else:
                    out_file.write("0\n")

snp_keep_file.close()
hap_file.close()
out_file.close()
