rm(list=ls())
library(tidyverse)
library(vegan)
library(reshape2)
library(pheatmap)

#### set up paths to input files
wd <- 'path/to/cwd'
meta_csv <- './data/results/metagenomics/metadata/sample_metadata_killuragh_3-6-22.forR.tsv'
taxid_dat <- './data/results/metagenomics/metadata/taxid_to_name.csv'
species_tab <- './data/results/metagenomics/OTU_tables/killuragh_teeth.dental_data_with_sources.25-04-22.species.tsv'
meta_csv <- '/home/iseult/Desktop/Killuragh/metatax_may_23/sample_metadata_killuragh_3-6-22.forR.tsv'
taxid_dat <- '/home/iseult/Desktop/Killuragh/metagenomics/updated_metatax_dataset_25-4-22/taxid_to_name.csv'
species_tab <- '/home/iseult/Desktop/Killuragh/metatax_may_23/killuragh_teeth.dental_data_with_sources.25-04-22.species.tsv'
# set working dir if necessary
#setwd(wd)

meta <- read.csv(meta_csv, sep='\t')
taxnames <- read.csv(taxid_dat,header=F)
colnames(taxnames) <- c("TAXID","Scientific_name","V3","V4","V5")
# remove unnecessary columns
taxnames['V3'] <- NULL
taxnames['V4'] <- NULL
taxnames['V5'] <- NULL

## Species level data
species_kraken <- read.csv(species_tab, sep='\t')
## convert raw counts to relative abundance
rownames(species_kraken) <- species_kraken$TAXID
species_kraken <- species_kraken[,-1]

MS <- as.matrix(species_kraken)

# Remove any columns of matrix with colSum < 1,000,000

MS_filt <- MS[,colSums(MS)>1000000]

MS1 <- sweep(MS_filt,2,colSums(MS_filt),"/")

# list of samples with colSums < 1 million reads colSums(MS) %>% sort()

## Some of these are really low abundance and not super informative 
# If proportion is less than 0.03% (so 0.00003?) throw it out

MS1[MS1 < 0.00003] <- 0
#get rid of rows with no counts left after filtering
MS1 <- MS1[rowSums(MS1[,-1])>0,]

## change taxonomy IDs to names
# convert normalised matrix back to dataframe

norm_sp <- data.frame(TAXID = row.names(MS1), MS1)
# convert taxon IDs to tax names

norm_sp_names <- merge(norm_sp,taxnames)
norm_sp_names['TAXID'] <- NULL

melt_sp <- melt(norm_sp_names)
colnames(melt_sp) <- c("Species", "ID","Proportion")  

# pathobionts
path_all <- melt_sp %>% filter(grepl("Streptococcus_mutans|Tannerella_forsythia|Treponema_denticola|Porphyromonas_gingivalis|Fusobacterium_nucleatum", Species))

# create matrix
pathM <- acast(path_all, ID~Species, value.var="Proportion")

# write this to a csv so you have a table
write.csv(pathM, "~/Desktop/relative_abundance_pathobionts.1M_read_threshold.csv")
# normalise matrix by max value in each column # general format: sweep(MS,2,colSums(MS),"/")
pathM_norm <- sweep(pathM,2,apply(pathM,2,max),"/")

row.order <- c("KGH1.A","KGH1.E","KGH2.B","KGH2.F","KGH6.A","ANN2.B","KK1","PN107","PN112","PN113","SRA62","CE003","KD026","KD045","KD046","KD057","X1H04","X1H06","X1H07","X1H13","X1H14","X2H10","X2H11","X2H17","LM_213_T","LM_306_T","LM_308_T","LM_309_T","LM_403_T","LM_406_T","LM_913_T","MDV248","HI1","HI2","HS1","HS2","HS3","L","O1","O4","Ajv70","Gok2","Gok5","Gok7","Gnie1","Kow41","Kow52","Kow58","Leg2","Leg6","Mas11","Mas15","Mas17","Mas19","Mas20","Mas22","Mas27","Mas7","Mas8","Niem28","Niem5","Niem7","Sow18","Sow8","Abusir1547t","Abusir1564t","Abusir1595t","Abusir1615t","Abusir1627t","Abusir1665t","Abusir1607t","Abusir1616t","Abusir1655t","Abusir1671t","Abusir1519c","Abusir1594c","DRT001","RUV001","SMD017","SMD046","SMD051","WIG001","SA.001","SA.002","SA.003","SA.004","SA.005","SA.006","SA.007","SA.008","SA.009","SA.010","SA.011","KT05","KT13","KT26","KT28","KT32","CS01","CS02","CS03","CS04","CS07","CS08","CS09","CS10","CS11","CS12","CS13","CS14","CS15","CS16","CS17","CS18","CS19","CS20","CS21","CS22","CS23","CS24","CS26","CS27","CS28","CS29","CS30","CS31","CS32","CS33","CS34","CS35","CS36","CS37","CS38","CS39","CS40","CS42","CS43","CS44","CS45","CS46","CS47","CS48","B61","G12","SRR061294","SRR061302","SRR062298","SRR062302","SRR346694","SRR514329","SRR059898","SRR059899","SRR061562","SRR061566","SRR062330","SRR353624","ERR1019366","ERR1022687","ERR1022692","ERR1034454","ERR1035437","ERR1035438","SRR059812","SRR061145","SRR061153","SRR061164","SRR061170","Ext_Ctls_IJ_NovaScreen","SRR1622808","SRR1622829","SRR1622855","SRR1622920","SRR1622921")

pathM_norm <- pathM_norm[row.order,]
pathM_norm <- pathM_norm[,c("Fusobacterium_nucleatum","Porphyromonas_gingivalis","Tannerella_forsythia","Treponema_denticola","Streptococcus_mutans")]
pheatmap(pathM_norm, fontsize_row=5, cellheight=3, cluster_cols=FALSE, cluster_rows=FALSE)

pathM <- pathM[row.order,]
pathM <- pathM[,c("Fusobacterium_nucleatum","Porphyromonas_gingivalis","Tannerella_forsythia","Treponema_denticola","Streptococcus_mutans")]

pheatmap(pathM, fontsize_row=5, cellheight=3, cluster_cols=FALSE, 
         cluster_rows=FALSE, display_numbers=TRUE, fontsize=3, fontsize_col=8)

