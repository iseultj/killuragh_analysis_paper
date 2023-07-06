#!/usr/bin/bash
set -eou pipefail
## note - need to set variables. Have put dummy variables as placeholders

# dependencies: megahit; scripts for concoct preparation; bwa; samtools


path_to_fastqs='path/to/fastqs'
list_of_fastqs='comma_separated_list'
mem=20000000000 # mem in b
threads=5
cut_up_fasta='/path/to/miniconda3/pkgs/concoct-1.1.0-py38h7be5676_2/bin/cut_up_fasta.py' # concoct script : concoct env
desman_lengths='/path/to/miniconda3/envs/desman_env/bin/Lengths.py' # script from desman
collate_covs='/path/to/miniconda3/envs/desman_env/bin/Collate.pl'
#a. MEGAHIT assembly
megahit -r ${list_of_fastqs}  --memory ${mem} -t ${threads} --out-dir megahit --out-prefix killuragh
wait
#b. CONCOCT preparation

#i) cut up to 10k length contigs (using concoct script)
python $cut_up_fasta -c 10000 -o 0 -m megahit/killuragh.contigs.fa > Assembly/killuragh.contigs_c10k.fa
wait
bwa index Assembly/killuragh.contigs_c10k.fa

#ii) map trimmed reads to contigs
for fq in $(ls ${path_to_fastqs}*gz); do id=$(echo $fq | rev | cut -d/ -f1 | rev | cut -d. -f1); bwa aln -l 16500 -n 0.01 -o 2 -t 10 Assembly/killuragh.contigs_c10k.fa ${fq} > mapping/${id}.sai 2> mapping/${id}.bwa.log; bwa samse Assembly/killuragh.contigs_c10k.fa mapping/${id}.sai ${fq} 2> mapping/${id}.samse.log | samtools view -Sb -F4 - 2> mapping/${id}.viewF4.log | samtools sort - -o mapping/${id}.sorted.bam 2> mapping/${id}.sort.log; done
wait

# iii) Calculate assembly contig lengths 
python $desman_lengths -i Assembly/killuragh.contigs_c10k.fa > Assembly/Lengths.txt
wait
# iv) Calculate coverages
for file in $(ls mapping/*bam); do id=$(echo $file | cut -d- -f1-2); bedtools genomecov -ibam ${file} -g Assembly/Lengths.txt > mapping/${id}_cov.txt; done
wait
for i in $(ls mapping/*_cov.txt); do echo $i; stub=${i%_cov.txt};  stub=${stub#Map\/}; echo $stub; awk -F"\t" '{l[$1]=l[$1]+($2 *$3);r[$1]=$4} END {for (i in l){print i","(l[i]/r[i])}}' $i > mapping/${stub}_cov.csv; done
wait
$collate_covs mapping > Coverage.csv

echo "Finished preparation for concoct: continue to next step!"
