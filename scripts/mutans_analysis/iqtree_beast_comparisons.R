rm(list=ls())


#libraries
library(tidyverse)
library(treeio)
library(ggtree)
library(data.table)
library(ggrepel)

wd <- './data/results/s_mutans/tree_data/'
beast_clade_df <- './beast_clades_median_mcc_tree.csv'
beast_dat <- './median_mcc_tree_beast_data.csv'
iqtree_treefile <-  '../iqtree/S_mutans.beast_sites.1000boots.treefile'

# 1. Get BEAST clade assignments


beast <- fread(beast_clade_df, data.table=F)
# 2. IQtree with names corrected to be same as beast

iqtree <- read.iqtree(iqtree_treefile)

# 3. Plot ML tree
ggtree(iqtree) + geom_nodelab(geom="label", node="internal") + theme_tree() # too messy to read
selected_nodes <- c(413,414,415,416,663,417,653,564,565,578,489,540,542,550,440, 664,775,666,756,776)
# just label high quality branches?
hiqual <- iqtree %>% as_tibble() %>% filter(UFboot > 95) %>% select(node)
tmp <- iqtree %>% as_tibble() %>% filter(UFboot > 95) %>% filter(branch.length > 0.002) 
#ggtree(iqtree)+ geom_text(aes(label=node),size=2) # get list of nodes you want to label

ggtree(iqtree) + geom_nodelab(data=td_filter(node %in% hiqual$node), geom="label", node="internal", size=2, 
                              label.padding = unit(0.2, "lines"),
                              label.r = unit(0.1, "lines"),label.size=0.2) + theme_tree() 

# 4. Plot ML tree with beast clade assignment as metadata matrix

rownames(beast) <- beast$ID
beast$ID <- NULL
b <- as.matrix(beast)
p <- ggtree(iqtree) + theme_tree()

p1 <- ggtree(iqtree) + geom_nodelab(data=td_filter(node %in% hiqual$node), geom="label", node="internal", size=2, 
                                    label.padding = unit(0.2, "lines"),
                                    label.r = unit(0.1, "lines"),label.size=0.2) + theme_tree() 

pdf("iqtree_with_beast_clades.confident_nodes_labelled.pdf", width = 10, height = 20)
gheatmap(p1, b, colnames=FALSE) + scale_x_ggtree()  +
  scale_fill_manual(values=c( "#7F3C8D", "#11A579",	"#3969AC",	"#F2B701",	
                                    "#E73F74","#BE0032", "#80BA5A","#A5AA99","black")) 
dev.off()
