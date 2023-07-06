#!/bin/bash

set -eou pipefail
### NOTE: you will need to change paths to programmes/reference genomes used in this script: they have been saved as variables to make this easier

bamlist=$1
path_to_gatk=$(echo /Software/GenomeAnalysisTK.jar)
path_to_reference=$(echo /Reference_Genomes/Pathogen_Genomes/T_forsythia/T_forsythia_92A2.NC_016610.fasta)
path_to_bcftools=$(echo /Software/bcftools-1.10.2/bcftools)
path_to_plink=$(echo plink) # note use plink1.9

# this is not parallelised for now

# Step 1: Unified Genotyper
# Note that --ploidy 2 was used for compatability with multivcfanalyzer, which was later excluded from analysis in the paper. 

while read line; do id=$(basename $line | cut -d. -f1); java -Xmx8g -jar ${path_to_gatk} -T UnifiedGenotyper -I ${line} -R ${path_to_reference} -o unified_genotyper_all/${id}.mq25_34bp.unifiedgt.vcf --output_mode EMIT_ALL_SITES -mbq 30 -ploidy 2 -gt_mode DISCOVERY; done < $bamlist

wait
# Step 2: merge vcfs

for file in $(ls unified_genotyper_all/*.mq25_34bp.unifiedgt.vcf); do bgzip $file; tabix -p vcf ${file}.gz; done

wait

ls unified_genotyper_all/*vcf.gz > vcfs_to_merge.txt

wait

${path_to_bcftools} merge --file-list vcfs_to_merge.txt -Oz -o T_forsythia_dataset.mq25.bp34.bq30.vcf.gz

wait

## Step 3: coverage filter

# Hard minimum 2X or 2 SDs below mean coverage; maximum 2SDs above max coverage. 

# this command only works for VCF in specific order (see data repository) - need to edit for further analysis.

$path_to_bcftools view T_forsythia_dataset.mq25.bp34.bq30.vcf.gz -e 'FORMAT/DP[0] > 158 & FORMAT/DP[0] < 13 & FORMAT/DP[1] < 17 & FORMAT/DP[1] > 155 & FORMAT/DP[2] < 15 & FORMAT/DP[2] > 155 & FORMAT/DP[3] < 15 & FORMAT/DP[3] > 155 &FORMAT/DP[4] < 15 & FORMAT/DP[4] > 155 & FORMAT/DP[5] < 11 & FORMAT/DP[5] > 157 & FORMAT/DP[6] < 11 & FORMAT/DP[6] > 157 & FORMAT/DP[7] < 11 & FORMAT/DP[7] > 157 & FORMAT/DP[8] < 11 & FORMAT/DP[8] > 156 & FORMAT/DP[9] < 2 & FORMAT/DP[9] > 99 & FORMAT/DP[10] < 2 & FORMAT/DP[10] > 47 & FORMAT/DP[11] < 2 & FORMAT/DP[11] > 19 & FORMAT/DP[12] < 2 & FORMAT/DP[12] > 28 & FORMAT/DP[13] < 2 & FORMAT/DP[13] > 16 & FORMAT/DP[14] < 2 & FORMAT/DP[14] > 24 & FORMAT/DP[15] < 2 & FORMAT/DP[15] > 14 & FORMAT/DP[16] < 2 & FORMAT/DP[16] > 13 & FORMAT/DP[17] < 2 & FORMAT/DP[17] > 12 & FORMAT/DP[18] < 2 & FORMAT/DP[18] > 9 &  FORMAT/DP[19] < 2 & FORMAT/DP[19] > 9 & FORMAT/DP[20] < 2 & FORMAT/DP[20] > 10 & FORMAT/DP[21] < 2 & FORMAT/DP[21] > 10 & FORMAT/DP[22] < 2 & FORMAT/DP[22] > 8 & FORMAT/DP[23] < 2 & FORMAT/DP[23] > 9 & FORMAT/DP[24] < 2 &FORMAT/DP[24] > 7' -Oz -o T_forsythia_dataset.mq25.bp34.bq30.covfilt.hardmin2.vcf.gz

wait

## Remove LowQual sites

${path_to_bcftools} view T_forsythia_dataset.mq25.bp34.bq30.covfilt.hardmin2.vcf.gz -f.,PASS -Oz -o T_forsythia_dataset.mq25.bp34.bq30.covfilt.hardmin2.rmlowqual.vcf.gz

wait 

## Remove sites with 100% missingness

$path_to_plink --vcf T_forsythia_dataset.mq25.bp34.bq30.covfilt.hardmin2.rmlowqual.vcf.gz --allow-extra-chr --geno 0.999 --make-bed --out T_forsythia_dataset.mq25.bp34.bq30.covfilt.hardmin2.rmlowqual.missfilt --const-fid

wait 
# remove invariant sites
$path_to_plink --bfile T_forsythia_dataset.mq25.bp34.bq30.covfilt.hardmin2.rmlowqual.missfilt --allow-extra-chr --maf 0.01 --allow-extra-chr --make-bed --out T_forsythia_dataset.mq25.bp34.bq30.covfilt.hardmin2.rmlowqual.missfilt.rm_invariant

wait

# prepare input for HapSNP script

cat T_forsythia_dataset.mq25.bp34.bq30.covfilt.hardmin2.rmlowqual.missfilt.rm_invariant.bim | awk '{print $1,$4-1,$4}' > T_forsythia.ascertained_min2x.pileup.bed

wait

cat T_forsythia_dataset.mq25.bp34.bq30.covfilt.hardmin2.rmlowqual.missfilt.rm_invariant.bim | awk '{print $1,$2,$3,$4,$5,$6}' > T_forsythia.ascertained_min2x.pileup.bim
wait

echo "ascertainment finished, don't forget to clean up intermediate files"


