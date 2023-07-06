#!/usr/bin/python
from __future__ import division
import sys
import os
import subprocess
from subprocess import call
import errno
from datetime import datetime
import csv
import argparse
from joblib import Parallel, delayed
import pysam
from multiprocessing import Queue

## use parser module
parser = argparse.ArgumentParser()
## add arguments to accept
parser.add_argument("-R","--ref",help="Reference genome")
parser.add_argument("-dcov","--depth_covg",help="filter minimum depth of coverage across all samples [0]",type=int,default=0)
parser.add_argument("-mem","--java_mem",help="Java memory (in GB) [25]",type=int,default=25)
parser.add_argument("-I","--input",help="list of bams")
parser.add_argument("-o","--output",help="output name",default="bacterial_heteroplasmy")
parser.add_argument("-s","--sites",help="variant sites in BED format to calculate contamination for")

## Collect Arguments
args = parser.parse_args()

path_to_genome = args.ref
dcov = str(args.depth_covg)
mem = str(args.java_mem)

bamlist  =  open(args.input, 'r')
all_files = ['.'.join(line.rstrip().split('.')[:-1]) for line in bamlist]
bamlist.close()
num_samples = len(all_files)

out=str(args.output)
final_file = open(out + ".contamination_estimates.txt",'w')
sites_file = str(args.sites)

path_to_gatk = "/Software/GenomeAnalysisTK.jar"
path_to_bcftools = "/Software/bcftools-1.10.2/bcftools"


for bam in all_files:
    i = bam.split('/')[-1]
    try :
        test = open(i + ".mutations.pileup",'r')
        print "Pileup already performed. Moving on..."
    except IOError:
        call("java -Xmx" + mem + "G -jar " + path_to_gatk  + " -T Pileup -R " + path_to_genome + " -I " + bam + ".bam -o " + i + ".mutations.pileup -L " + sites_file, shell=True)  
# Filter for Base quality 30
    try :
        test = open(i + ".mutations.bq30.pileup",'r')
        print "Base quality filtering already performed. Moving on..."
    except IOError:
        input_file_name= i + ".mutations.pileup"
        output_file_name= i + ".mutations.bq30.pileup"
        input_file= open(input_file_name,'r')
    #for input pileup file, ignore header and footer lines
        for line in input_file:
            if  "REDUCE" in line:
                pass

    # Filter pileup file by base quality 30
            else:
                each = line.strip()
                each = each.split(' ')
                my_SNP = each[3]
                quality = each[4]
                quality = quality.strip()

                counter = 0
                quality_list = []
                mySNP_list = []

    # Cycle through each value in quality string. Quality offset is 33.
                for q in quality[0:len(quality)]:
                    offset_quality = ord(q)-33

    # If the quality of the base is greater than or equal to 20, store this base and the corresponding quality in a list (well 2 lists - one for bases, one for qualities.)
                    if offset_quality >= 30:
                        quality_list.append(q)
                        mySNP_list.append(my_SNP[counter])
            
                    counter = counter + 1

    # Output high quality SNP calls to a file. 	
                if len(quality_list) > 0:
                    header_string = '\t'.join(each[0:3])
                    mySNP_string = ''.join(mySNP_list)
                    quality_string = ''.join(quality_list)
                    output_filtered_file = open(output_file_name,'a')
                    output_filtered_file.write(header_string + '\t' + mySNP_string + '\t' + quality_string + '\n')
        output_filtered_file.close()

    #Raw Base Counting
    try : 
        test = open(i + ".mutations.bq30.basecounts.txt",'r')
        print "Base Counting already performed. Moving on..."
    except IOError:
        base_count_input= open(i + ".mutations.bq30.pileup",'r')
        base_count_output = open(i + ".mutations.bq30.basecounts.txt",'w')
        base_count_output.write("Position Ref_Allele A C G T" + '\n')
        for line in base_count_input :
            each = line.strip()
            each = each.split()
            calls=each[3]
            position=each[1]
            ref_allele = each[2]
            A = calls.count('A')
            C = calls.count('C')			
            G = calls.count('G')
            T = calls.count('T')
            base_count_output.write(position + " " + ref_allele  + " " + str(A)  + " " + str(C)  + " " + str(G)  + " " + str(T) + '\n')
        base_count_output.close()

    #Major minor allele identification

    try:
        test = open(i + ".mutations.bq30.major_minor.txt",'r')
        print "Major/minor allele identification already performed. Moving on..."

    except IOError:
        major_minor_input= open(i + ".mutations.bq30.basecounts.txt",'r')
        major_minor_output = open(i + ".mutations.bq30.major_minor.txt",'w')
        major_minor_output.write("Position Ref_Allele Major_Allele Major_Allele_Count Minor_Alleles_Count Minor_Alleles_Count_no_MD" + '\n')
        total_major_alleles = []
        total_major_alleles_hetplas_sites_only =[]
        total_minor_alleles = []
        total_minor_alleles_hetplas_sites_only = []
        total_minor_alleles_noMD = []
        total_minor_alleles_noMD_hetplas_sites_only = []	
        for line in major_minor_input :
            if "Ref_Allele" in line:
                pass
            else:
                each = line.strip()
                each = each.split()
                position = each[0]
                ref_allele = each[1]
                alleles = map(int, each[2:6])
                major_allele_count = max(alleles)
                if  alleles.index(major_allele_count) is 0:
                    major_allele = "A"
                    minor_allele_count = sum(alleles[1:4])
                    minor_allele_count_noMD = sum(alleles[1:4])
                if  alleles.index(major_allele_count) is 1:
                    major_allele = "C"
                    minor_allele_count = sum(alleles[0:1] + alleles[2:4])
                    minor_allele_count_noMD = sum(alleles[0:1] + alleles[2:3])
                if  alleles.index(major_allele_count) is 2:
                    major_allele = "G"
                    minor_allele_count = sum(alleles[0:2] + alleles[3:4])
                    minor_allele_count_noMD = sum(alleles[1:2] + alleles[3:4])
                if  alleles.index(major_allele_count) is 3:
                    major_allele = "T"
                    minor_allele_count = sum(alleles[0:3])
                    minor_allele_count_noMD = sum(alleles[0:3])
                major_minor_output.write(position + " " + ref_allele + " " + major_allele + " " + str(major_allele_count) + " " + str(minor_allele_count) + " " + str(minor_allele_count_noMD) + '\n')
                total_major_alleles.append(major_allele_count)
                total_minor_alleles.append(minor_allele_count)
                total_minor_alleles_noMD.append(minor_allele_count_noMD)
                if minor_allele_count > 0 :
                    total_major_alleles_hetplas_sites_only.append(major_allele_count)
                    total_minor_alleles_hetplas_sites_only.append(minor_allele_count)
                    total_minor_alleles_noMD_hetplas_sites_only.append(minor_allele_count_noMD)					
        major_minor_output.close()
	try:
		all_sites_contam = str(100*(sum(total_minor_alleles)/(sum(total_major_alleles)+sum(total_minor_alleles))))
	except ZeroDivisionError:
		all_sites_contam = '0'

	try:
		all_sites_no_md_contam = str(100*(sum(total_minor_alleles_noMD)/(sum(total_major_alleles)+sum(total_minor_alleles_noMD))))
	except ZeroDivisionError:
		all_sites_no_md_contam = '0'
	try:
		hetplas_sites_contam = str(100*(sum(total_minor_alleles_hetplas_sites_only)/(sum(total_major_alleles_hetplas_sites_only)+sum(total_minor_alleles_hetplas_sites_only))))
	except ZeroDivisionError:
		hetplas_sites_contam = '0'
	try:
		hetplas_sites_nomd_contam = str(100*(sum(total_minor_alleles_noMD_hetplas_sites_only)/(sum(total_major_alleles_hetplas_sites_only)+sum(total_minor_alleles_noMD_hetplas_sites_only))))
	except ZeroDivisionError:
		hetplas_sites_nomd_contam = '0'

        final_file.write(i + " " + all_sites_contam + " " + all_sites_no_md_contam + " " + hetplas_sites_contam + " " + hetplas_sites_nomd_contam + '\n') 
final_file.close()
