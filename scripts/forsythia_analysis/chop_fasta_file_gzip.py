#!/usr/bin/python

# Quick script to read in fasta file (saving header)
# Split sequence to kmers with n-bp slide (these are args)
#Turn this into a dummy fastq to align to coding seqs
# Output is printed to stdout: pipe to bgzip and redirect to fastq
# From this align to gene ref seqs.

from __future__ import division
import argparse
import os
import sys
import gzip

parser = argparse.ArgumentParser()

parser.add_argument("-f","--fasta",help="input fasta file to split")
parser.add_argument("-k","--k",help="length of kmer to split fasta file to")
parser.add_argument("-n","--n",help="Length in bp of slide")

args = parser.parse_args()

count = 0
f = gzip.open(args.fasta,"rt")
name = str(args.fasta).split(".")[0]
k = int(args.k)
n = int(args.n)
m = k-n
dummy_qual = "M"*k
line = f.readline()
#header 1
header = "@" + line.strip()[1:]
print("first header",file=sys.stderr,flush=True)
#rest of file
line = f.readline()
#string to save contig seq to 
s = ''
#count
count=0
linecount=0
print("Starting processing: Split to reads of length " + str(k) + " slide of length " + str(n) + "." + "\n", file=sys.stderr,flush=True)
while line:
	if line.strip().startswith(("A","C","T","G","N","a","c","t","g","n")):
		s = s + str(line.strip())
		linecount+=1
		print(str(linecount) + " sequence lines recorded", file=sys.stderr,flush=True)
		line = f.readline()
	elif line.strip().startswith(">"):
		#new contig
		#do processing first
		for i in range(len(s)-m):
			seq = s[i:i+k]
			if seq=="N"*k:
				print("skipped " + str(i) + "th subseq",file=sys.stderr)
				continue
			elif seq=="n"*k:
				print("skipped " + str(i) + "th subseq",file=sys.stderr)
				continue
			else:
				print(header + "_" + str(i))
				print(str(seq))
				print("+")
				print(dummy_qual)
				count +=1

		#now move on to next read
		header = "@" + line.strip()[1:]
		s = ''
		line = f.readline()
#final contig

for i in range(len(s)-m):
	seq = s[i:i+k]
	if seq=="N"*k:
		print("skipped " + str(i) + "th subseq",file=sys.stderr,flush=True)
		continue
	elif seq=="n"*k:
		print("skipped " + str(i) + "th subseq",file=sys.stderr,flush=True)
		continue
	else:
		print(header + "_" + str(i))
		print(str(seq))
		print("+")
		print(dummy_qual)
		count +=1

f.close()
print(str(count) + " dummy reads written!",file=sys.stderr)
