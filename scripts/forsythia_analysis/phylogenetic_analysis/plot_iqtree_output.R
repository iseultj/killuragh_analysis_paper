rm(list=ls())

wd <- './data/results/forsythia/phylogenetic_analysis/iqtree_neanderthal_outgroup/'
# plot ML trees forsythia

setwd(wd)

library(tidyverse)
library(ggtree)
library(treeio)
library(ggrepel)
library(ggbreak)

# plot independent runs to check same/similar

tr <- read.iqtree('run1/T_forsythia.ind40.keep_special.geno5.seed_3495.B1000.iqtree.timetree.nwk')
tr <- read.iqtree('run2/T_forsythia.ind40.keep_special.geno5.seed_31065.B1000.iqtree.timetree.nwk')
tr <- read.iqtree('run3/T_forsythia.ind40.keep_special.geno5.seed_8614.B1000.iqtree.timetree.nwk')
tr <- read.iqtree('run4/T_forsythia.ind40.keep_special.geno5.seed_9099.B1000.iqtree.timetree.nwk')


ggtree(tr) + geom_tiplab() + geom_nodelab(geom="label") + theme_tree2() # + xlim(0,0.6)


# # Dog outgroup

tr_dog <- read.iqtree('../iqtree_dog_outgroup/run1/T_forsythia.2x_ascertain.strict_hapsnp.clip_nonudg.bq30.mq25.34bp.mind40.keep_special.seed_23246.B1000.iqtree.timetree.nwk')
tr_dog <- read.iqtree('../iqtree_dog_outgroup/run2/T_forsythia.2x_ascertain.strict_hapsnp.clip_nonudg.bq30.mq25.34bp.mind40.keep_special.seed_32008.B1000.iqtree.timetree.nwk')
tr_dog <- read.iqtree('../iqtree_dog_outgroup/run3/T_forsythia.2x_ascertain.strict_hapsnp.clip_nonudg.bq30.mq25.34bp.mind40.keep_special.seed_4478.B1000.iqtree.timetree.nwk')
tr_dog <- read.iqtree('../iqtree_dog_outgroup/run4/T_forsythia.2x_ascertain.strict_hapsnp.clip_nonudg.bq30.mq25.34bp.mind40.keep_special.seed_13683.B1000.iqtree.timetree.nwk')

ggtree(tr_dog) + geom_tiplab() + geom_nodelab(geom="label") + theme_tree2() #+ xlim(0,0.4)
