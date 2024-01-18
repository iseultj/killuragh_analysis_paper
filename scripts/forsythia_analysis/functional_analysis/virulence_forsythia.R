rm(list=ls())

wd <-- './data/results/forsythia/functional_analysis/'
setwd(wd)

library(tidyverse)
library(gplots)
library(data.table)
df <- fread('virulence_factors_all.txt', data.table=F)

colnames(df) <- c("chr","start","end","reads_overlapping_interval",
"bases_overlapping_interval","length_of_interval","fraction_covered","id","tmp")
sample_meta <- fread('sample_nreads.txt', data.table=F)
df <- merge(df, sample_meta)
gene_meta <- fread('gene_meta.txt', data.table=F)
df <- merge(df, gene_meta)
#df <- df[order(df$category),]
df$Gene <- factor(df$Gene, levels=c("nanH","siaHI","susB","hexA","mgsA","oxyR","dppIV","BFO_RS01890",
                                    "prtH","miropin","pse","tfsA","tfsB","porU","porT","sov","porK",
                                    "bspA","bspB","karilysin","miropsin-1","miropsin-2","mirolase","mirolysin","forsilysin"))
# filter out bad samples
df <- df %>% filter(id !="OME003") %>% filter(id != "CS07") %>% filter(id != "CS09") %>% filter(id != "CS10") %>% filter(id != "CS12") %>% filter(id != "CS39")


# normalise reads
df$normalised_reads <- df$reads_overlapping_interval/(df$nreads_total*(df$length_of_interval/3405521))

#decide if present or absent based on normalised reads > 0.8, more than 50% covered

df <- df %>% mutate(presence = if_else(fraction_covered > 0.95, 1,if_else((normalised_reads > 0.8 & fraction_covered > 0.5),1,0)))

ggplot(df, aes(x=id, y=Gene, fill=presence)) + geom_tile() +
  scale_fill_gradient2(low="white", midpoint=0.5, mid="#E2D8D3FF", high="black") +
  theme_classic() + theme(axis.text.x = element_text(angle=90)) +ggtitle("Virulence factor presence/absence rough plot")


# sanity check - read counts vs coverage?

s <- df %>% select(id,presence)
counts <- s %>% group_by(id)  %>% table() %>% as.data.frame() %>% filter(presence == 1)
cov <- fread('coverage_data.tab', data.table=F)
colnames(cov)
colnames(counts) <- c("ID","Presence","Count")
counts <- merge(counts,cov)
ggplot(counts, aes(x=Coverage,y=Count)) + geom_point() + theme_bw() # facet by cov -anything above 22 is modern

counts <- counts %>% mutate(date=if_else(Coverage > 22,"Modern",if_else(ID=="VLC009","Modern",if_else(ID=="VLC004","Modern","Ancient"))))

library(ggrepel)
ggplot(counts, aes(x=Coverage,y=Count)) + geom_point() + theme_bw() + facet_grid(~date)  + geom_text_repel(aes(label=ID)) + ggtitle("Sanity check - T. forsythia coverage vs virulence counts")
                  


# keep in order of tree



df$id <- factor(df$id, levels=(c("OH2617_COT023","GOY005","KGH1","KGH2","2H10","TAF008","OAK005","G12",
                                 "SMD046","SMD051","WIG001","LM_213_T","CS23","CS37",
                                 "CS31","CS32","UB22","VLC009","ATCC43037",
                                 "WW10960","3313","KS16","9610","UB4","WW11663","UB20","VLC004")))


ggplot(df, aes(y=id, x=Gene, fill=presence)) + geom_tile() +
  scale_fill_gradient2(low="white", midpoint=0.5, mid="#E2D8D3FF", high="black") +
  theme_classic() + theme(axis.text.x = element_text(angle=90)) +ggtitle("Virulence factor presence/absence")

# Exclude dog for final plot
df <- df %>% filter(ID != "OH2617_COT023")

df$id <- factor(df$id, levels=(c("GOY005","KGH1","KGH2","2H10","TAF008","OAK005","G12",
                                 "SMD046","SMD051","WIG001","LM_213_T","CS23","CS37",
                                 "CS31","CS32","UB22","VLC009","ATCC43037",
                                 "WW10960","3313","KS16","9610","UB4","WW11663","UB20","VLC004")))

# Reorder genes for clarity in figure
df$Gene <- factor(df$Gene, levels=c("oxyR","nanH","siaHI","susB","hexA","mgsA","dppIV","BFO_RS01890",
                                    "pse","tfsA","tfsB","porU","porT","sov","porK",
                                    "bspA","bspB","prtH","karilysin","miropsin-1","miropsin-2",
                                    "mirolase","mirolysin","forsilysin","miropin"))

ggplot(df, aes(y=id, x=Gene, fill=presence)) + geom_tile() +
  scale_fill_gradient2(low="white", midpoint=0.5, mid="#E2D8D3FF", high="black") +
  theme_classic() + theme(axis.text.x = element_text(angle=90)) +ggtitle("Virulence factor presence/absence")

