#!/bin/python
#convert vcf to fasta
# Conor Rossi 2022 with modifications by Iseult Jackson
#1. vcf-to-tab < merged.vcf > merged_tabbed.vcf

#2. then run this script i.e. python vcf2fasta.py


import sys, os, io
n_samples=int(sys.argv[1])

with open('merged_tabbed.vcf') as f:
	for i in range(3,n_samples):
		first_line = f.readline()
		first_line = first_line.strip().split()
		first_line = first_line[i]
		first_line = first_line.split("_")[0]
		f.seek(0)
		fasta = open(first_line + ".fa",'w')
		fasta.write(">" + first_line + "\n")
		
		next = f.readlines()[1:]
		for lines in next:
			field = lines.strip().split()
			geno = field[i]
			base = geno[0]
			if base == ".":
				fasta.write("N")
			else:
				fasta.write(base)
		fasta.write("\n")
		f.seek(0)
