import os, sys

out_file=open(sys.argv[1], "w")
n=int(sys.argv[2])

for i in range(0,n):
    out_file.write("0 0 ")
for i in range(0, 99):
    out_file.write("1 1 ")
for i in range(0, 108):
    out_file.write("2 2 ")
out_file.write("\n")
out_file.close()
