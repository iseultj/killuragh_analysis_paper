rm(list=ls())

library(tidyverse)
library(data.table)
# 
wd <- './data/results/denticola/heteroplasmy/heteroplasmy_estimates_competitive/chromosome_level_treponemes/'
# need to change wd for different alignments
setwd(wd)
positions <- fread('major_minor_all_positions.txt', data.table = F)
# use data.frame() so they've different memory addresses
major_allele_count <- data.frame(positions)
minor_allele_count <- data.frame(positions)
minor_allele_count_no_md <- data.frame(positions)
indfiles <- list.files(pattern="*mutations.bq30.major_minor.txt.gz", include.dirs=TRUE)
for (i in indfiles) {
  tmp <- fread(i)
  tmp_major_allele_count <- tmp %>% select("Position","Ref_Allele","Major_Allele_Count")
  tmp_major_allele_count$Major_Allele_Count <- as.numeric(tmp_major_allele_count$Major_Allele_Count)
  colnames(tmp_major_allele_count) <- c("Position","Ref_Allele",i)
  tmp_minor_allele_count <- tmp %>% select("Position","Ref_Allele","Minor_Alleles_Count")
  tmp_minor_allele_count$Minor_Alleles_Count <- as.numeric(tmp_minor_allele_count$Minor_Alleles_Count)
  colnames(tmp_minor_allele_count) <- c("Position","Ref_Allele",i)
  tmp_minor_allele_count_no_md <- tmp %>% select("Position","Ref_Allele","Minor_Alleles_Count_no_MD")
  tmp_minor_allele_count_no_md$Minor_Alleles_Count_no_MD <- as.numeric(tmp_minor_allele_count_no_md$Minor_Alleles_Count_no_MD)
  colnames(tmp_minor_allele_count_no_md) <- c("Position","Ref_Allele",i)
  major_allele_count <- merge(major_allele_count, tmp_major_allele_count,by.x=c("Position","Ref_Allele"),by.y=c("Position","Ref_Allele"), all.x=TRUE)

  minor_allele_count <- merge(minor_allele_count, tmp_minor_allele_count,by.x=c("Position","Ref_Allele"),by.y=c("Position","Ref_Allele"), all.x=TRUE)

  minor_allele_count_no_md <- merge(minor_allele_count_no_md, tmp_minor_allele_count_no_md,by.x=c("Position","Ref_Allele"),by.y=c("Position","Ref_Allele"), all.x = TRUE)

}

# Aggregate all of this data


rownames(major_allele_count) <- major_allele_count$Position
ref <- major_allele_count %>% select("Position","Ref_Allele")
major_allele_count$Position <- NULL
major_allele_count$Ref_Allele <- NULL
# set NA to 0
major_allele_count[is.na(major_allele_count)] <- 0
major_sums <- rowSums(major_allele_count)
major_sums <- data.frame(major_sums)
major_sums$Position <- rownames(major_sums)
major_sums <- merge(major_sums, ref)
colnames(major_sums)
# minor
rownames(minor_allele_count) <- minor_allele_count$Position
ref <- minor_allele_count %>% select("Position","Ref_Allele")
minor_allele_count$Position <- NULL
minor_allele_count$Ref_Allele <- NULL
# set NA to 0
minor_allele_count[is.na(minor_allele_count)] <- 0
minor_sums <- rowSums(minor_allele_count)
minor_sums <- data.frame(minor_sums)
minor_sums$Position <- rownames(minor_sums)
minor_sums <- merge(minor_sums, ref)

# minor no md
rownames(minor_allele_count_no_md) <- minor_allele_count_no_md$Position
ref <- minor_allele_count_no_md %>% select("Position","Ref_Allele")
minor_allele_count_no_md$Position <- NULL
minor_allele_count_no_md$Ref_Allele <- NULL
# set NA to 0
minor_allele_count_no_md[is.na(minor_allele_count_no_md)] <- 0
minor_sums_no_md <- rowSums(minor_allele_count_no_md)
minor_sums_no_md <- data.frame(minor_sums_no_md)
minor_sums_no_md$Position <- rownames(minor_sums_no_md)
minor_sums_no_md <- merge(minor_sums_no_md, ref)

sums <- merge(major_sums, minor_sums)
sums <- merge(sums, minor_sums_no_md)

write.csv(sums, "aggregated_maj_minor_sums.comp_aln.csv")
