rm(list=ls())

wd <- './data/results/phylogenetic_analysis/beast_geno5_allmedieval_on/partitioned/'
setwd(wd)
library(treeio)
library(ape)
library(ggtree)
library(ggplot2)
library(data.table)
library(TreeTools)
library(cowplot)

beast1 <- read.beast('T_forsythia.medieval_on_no_dog.geno05.rm_invariant.partitioned.seed24081.burnin10pc.min_posterior90pc.mcc.nex')
beast2 <- read.beast('T_forsythia.medieval_on_no_dog.geno05.rm_invariant.partitioned.seed22913.burnin10pc.min_posterior90pc.mcc.nex')


t1 <- ggtree(beast1)  +
  geom_range(range='CAheight_0.95_HPD', color='#537d90', alpha=.5, size=1) +
  geom_nodelab(aes(x=branch, label=round(posterior, 2),geom="label"), vjust=-.5, size=3) + geom_tiplab(align=TRUE) +
  theme_tree2()

t2 <- ggtree(beast2)  +
  geom_range(range='CAheight_0.95_HPD', color='#537d90', alpha=.5, size=1) +
  geom_nodelab(aes(x=branch, label=round(posterior, 2),geom="label"), vjust=-.5, size=3) + geom_tiplab(align=TRUE) +
  theme_tree2()
revts(t1)
revts(t2)

bsky1 <- fread('seed_24081.bsky_plot.tab', data.table = F, skip = 1)
bsky2 <- fread('seed_22913.bsky_plot.tab', data.table = F, skip = 1)
head(bsky1)
ggplot(bsky1, aes(x=-time, y=log(median))) + geom_line() + geom_line(aes(x=-time, y=log(upper)), linetype="dotted") +
  geom_line(aes(x=-time,y=log(lower)), linetype="dotted") + theme_bw() + 
  scale_x_continuous(breaks = round(seq(-max(bsky1$time), -min(bsky1$time), by = 50),50))

ggplot(bsky2, aes(x=-time, y=log(median))) + geom_line() + geom_line(aes(x=-time, y=log(upper)), linetype="dotted") +
  geom_line(aes(x=-time,y=log(lower)), linetype="dotted") + theme_minimal()
