# Killuragh: relative abundance plots
rm(list=ls())
library(tidyverse)
library(vegan)
library(reshape2)

#### set up paths to input files
wd <- 'path/to/cwd'
meta_csv <- './data/results/metagenomics/metadata/sample_metadata_killuragh_3-6-22.forR.tsv'
taxid_dat <- './data/results/metagenomics/metadata/taxid_to_name.csv'
species_tab <- './data/results/metagenomics/OTU_tables/killuragh_teeth.dental_data_with_sources.25-04-22.species.tsv'
genus_tab <- './data/results/metagenomics/OTU_tables/killuragh_teeth.dental_data_with_sources.25-04-22.genus.tsv'

### Run script for plotting
# set working dir if necessary (uncomment line)
#setwd(wd)

meta <- read.csv(meta_csv, sep='\t')
taxnames <- read.csv(taxid_dat,header=F)
colnames(taxnames) <- c("TAXID","Scientific_name","V3","V4","V5")
# remove unnecessary columns
taxnames['V3'] <- NULL
taxnames['V4'] <- NULL
taxnames['V5'] <- NULL
#read in dataframes
species_kraken <- read.csv(species_tab, sep='\t')
genus_kraken <- read.csv(genus_tab, sep='\t')
#Clean up initial dfs
rownames(species_kraken) <- species_kraken$TAXID
species_kraken <- species_kraken[,-1]
rownames(genus_kraken) <- genus_kraken$TAXID
genus_kraken <- genus_kraken[,-1]
taxnames$TAXID <- as.character(taxnames$TAXID)
taxnames$Scientific_name <- as.character(taxnames$Scientific_name)
#Normalise counts. Just work with genus from here on out
G <- as.matrix(genus_kraken)
G <- sweep(G,2,colSums(G),"/")
# remove those with < 5%
G[G < 0.05] <- 0
#get rid of rows with no counts left after filtering
G <- G[rowSums(G[,-1])>0,]
norm_genus <-  data.frame(TAXID = row.names(G), G)
norm_genus_names <- merge(norm_genus,taxnames)
norm_genus_names['TAXID'] <- NULL

# melt 
gmelt <- melt(norm_genus_names)
colnames(gmelt) <- c("Scientific_name","ID","Rel_Ab")
# merge with metadata
gmelt_meta <- merge(gmelt,meta)
str(gmelt_meta)
# Drop all unecessary columns
gmelt_meta$Read_Depth <- NULL
gmelt_meta$Path_to_bracken_report <- NULL
gmelt_meta$SourceSink <- NULL
gmelt_meta$Publication <- NULL
gmelt_meta$Period <- NULL
gmelt_meta$Env <- NULL
gmelt_meta$Country <- NULL
# per-genus averages across bins: cast then re-melt for plotting
avgs_genus <- acast(gmelt_meta, Scientific_name ~ SampleSource, value.var = "Rel_Ab", mean)
# Add row for low abundance proportion
low_abundance <- 1-colSums(avgs_genus)
avgs_genus <- rbind(avgs_genus, low_abundance)

avg_genus_melt <- melt(avgs_genus)
colnames(avg_genus_melt) <- c("Scientific_Name","Source","Proportion")
str(avg_genus_melt)


pdf("killuragh_genus_plot.5pc_cutoff.pdf", width=40, height=20)

ggplot(avg_genus_melt, aes(x=Source, y=Proportion, 
                           fill=Scientific_Name, color=Scientific_Name)) +  
  geom_bar(position="fill",stat="identity") + theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

dev.off()


# Plot relative abundance of streptococcal species

# Streptococcal species **only**

# Dataframe: species_kraken
# Calculate relative abundance

S <- as.matrix(species_kraken)
S <- sweep(S,2,colSums(S),"/")

# For now, don't bother removing low-abundance- can remove them later.
norm_species <-  data.frame(TAXID = row.names(S), S)
norm_species_names <- merge(norm_species,taxnames)
norm_species_names['TAXID'] <- NULL

# filter for streprococci

strep <- norm_species_names %>% filter(grepl("Streptococcus", Scientific_name)) %>% filter(!grepl("phage", Scientific_name))

# convert to matrix and re-do relative abundance wrt Strep species

rownames(strep) <- strep$Scientific_name
strep <- strep[,1:396] # drop scientific name so matrix is numeric
strepM <- as.matrix(strep)
# remove any with relative abundance less than 0.5%
strepM[strepM < 0.005] <- 0
#get rid of rows with no counts left after filtering
strepM <- strepM[rowSums(strepM[,-1])>0,]
colSums(strepM)
# re-do rel ab
strepM1 <- sweep(strepM, 2, colSums(strepM),"/")
# Replace NA values with zeros
strepM1[is.na(strepM1)] <- 0

colSums(strepM1)
norm_strep <-  data.frame(Scientific_name = row.names(strepM1), strepM1)

# melt
smelt <- melt(norm_strep)
colnames(smelt) <- c("Scientific_name","ID","Rel_Ab")
# merge with metadata
smelt_meta <- merge(smelt,meta)
str(smelt_meta)
# Drop all unecessary columns
smelt_meta$Read_Depth <- NULL
smelt_meta$Path_to_bracken_report <- NULL
smelt_meta$SourceSink <- NULL
smelt_meta$Publication <- NULL
smelt_meta$Period <- NULL
smelt_meta$Env <- NULL
smelt_meta$Country <- NULL
# per-strep species averages across bins: cast then re-melt for plotting
avgs_strep <- acast(smelt_meta, Scientific_name ~ SampleSource, value.var = "Rel_Ab", mean)
#no strep

no_strep <- 1-colSums(avgs_strep)
avgs_strep <- rbind(avgs_strep, no_strep)

avgs_strep_melt <- melt(avgs_strep)
colnames(avgs_strep_melt) <- c("Scientific_Name","Source","Proportion")

ggplot(avgs_strep_melt, aes(x=Source, y=Proportion, fill=Scientific_Name, color=Scientific_Name)) +  
  geom_bar(position="fill",stat="identity") +
  theme_classic() + theme(axis.text.x=element_text(angle=90))

