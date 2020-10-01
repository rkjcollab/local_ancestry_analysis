import os, sys

out_file=open(sys.argv[1], "w")
n_admixed=int(sys.argv[2])
n_european=int(sys.argv[3])
n_african=int(sys.argv[4])

for i in range(0,n_admixed):
    out_file.write("0 0 ")
for i in range(0, n_european):
    out_file.write("1 1 ")
for i in range(0, n_african):
    out_file.write("2 2 ")
out_file.write("\n")
out_file.close()
