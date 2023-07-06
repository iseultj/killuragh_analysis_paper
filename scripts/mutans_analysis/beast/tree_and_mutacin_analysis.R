rm(list=ls())

library(treeio)
library(ape)
library(ggtree)
library(ggplot2)
library(data.table)
library(TreeTools)
library(cowplot)
library(tidyverse)

# set up paths:
wd <- './data/results/s_mutans/tree_data/'
beast_tree <- '../beast/S_mutans_aln.snp-sites.rm_multi.geno05.bmod.bcoalsky.combined.1e8_burnin.mcc_no_min_posterior.median_heights.nex'
mutacin_data <- './mutacin_presence_absence.csv'
beast_clade_df <- './beast_clades_median_mcc_tree.csv'

beast_dat <- './median_mcc_tree_beast_data.csv'

metadata_tab <- './full_metadata.corrected_plink.23-06-22.csv'
lat_long_df <- './lat_long_countries.csv'
setwd(wd)

beast <- read.beast(beast_tree)
beast_dat <- as.tibble(beast)
# save this if necessary 
#fwrite(beast_dat, "beast_data.mcc_median.no_min_posterior.csv")

# exploratory plots
t1 <-   ggtree(beast) + theme_tree2()
revts(t1)

ggtree(beast)  +
  geom_nodelab(aes(x=branch, label=round(posterior, 2),geom="label"), vjust=-.5, size=3) + 
  ggtitle("S. mutans biallelic SNP tree; 1e8 burnin")


# read in as phy object for annotation data.

tr <- read.nexus(beast_tree)

tr_order <- Preorder(tr)

# organise by clade
t1 <- Subtree(tr_order, 812)
t2 <- Subtree(tr_order, 791)
t3 <- Subtree(tr_order, 798)
t4 <- Subtree(tr_order, 688)
t5 <- Subtree(tr_order, 679)
t6 <- Subtree(tr_order, 601)
t7 <- Subtree(tr_order, 415)
t8 <- Subtree(tr_order, 598) # this has pretty poor support - 
#not a clade, just some long branches within 414 clade which group together

# Put into dataframe and save this info

beast_clades <- rbind(cbind(t1[["tip.label"]],rep("clade1",length(t1[["tip.label"]]))),
                      cbind(t2[["tip.label"]],rep("clade2",length(t2[["tip.label"]]))),
                      cbind(t3[["tip.label"]],rep("clade3",length(t3[["tip.label"]]))),
                      cbind(t4[["tip.label"]],rep("clade4",length(t4[["tip.label"]]))),
                      cbind(t5[["tip.label"]],rep("clade5",length(t5[["tip.label"]]))),
                      cbind(t6[["tip.label"]],rep("clade6",length(t6[["tip.label"]]))),
                      cbind(t7[["tip.label"]],rep("clade7",length(t7[["tip.label"]]))),
                      cbind(t8[["tip.label"]],rep("clade8",length(t8[["tip.label"]]))))

tr_grouped <- groupClade(tr_order, c(812,791,798,688,679,601,415,598))



options(ignore.negative.edge=F)
tree_coloured <- ggtree(tr_grouped, aes(color=group)) + theme(legend.position='none') + coord_cartesian(clip = 'off') + 
  theme_tree2(plot.margin=margin(6, 120, 6, 6)) +
  theme(legend.position=c(.1, .8))  +
  scale_color_manual(values=c("black", "#7F3C8D", "#11A579",	"#3969AC",	"#F2B701",	
                              "#E73F74","#BE0032", "#80BA5A","#A5AA99","#008695",		"#F97B72"))

# put killuragh on top
tree_coloured <- tree_coloured  %>% rotate(412) %>% rotate(686) %>% rotate(687)
revts(tree_coloured)
tip_order <- tree_coloured$data %>% select("label", "y") # tip order the same as when you have positive edge lengths. 
#fwrite(tip_order, "tip_order.mcc_tree_rotated.medians.set_negative_edge_to_zero.no_min_posterior.24-05-23.csv")
m <- fread(mutacin_data, data.table=F)

rownames(m) <- m$ID
m$ID <- NULL
M <- as.matrix(m)
# tip order from tree
tip_order <- c("Cluster10","GCF-000229945.1","GCF-000229285.1","AHRY01-N29","GCF-020531565.1",
               "GCF-012641575.1","GCF-000522745.1","AHRZ01-NMT4863","GCF-020529995.1","GCF-000522885.1",
               "AHSG01-NFSM1","GCF-012641605.1","ERS537508","SRS2403914","SRS841386","GCF-020529885.1",
               "GCF-020530795.1","GCF-020529845.1","SRS3070134","GCF-015556125.1","GCF-020530255.1",
               "GCF-020531165.1","GCF-020530965.1","GCF-020531185.1","GCF-000228785.1","GCF-012641765.1",
               "AHRP01-11SSST2","GCF-000229405.1","GCF-000228805.1","GCF-000230225.1","AHRU01-5SM3",
               "GCF-020530545.1","AHRD01-SA38","AHRL01-4SM1","GCF-000229445.1","GCF-000229485.1",
               "GCF-012642235.1","GCF-020531545.1","GCF-002083175.1","GCF-002083175.2","AHRK01-1SM1",
               "ERS537544","GCF-012641545.1","GCF-012642005.1","GCF-020531625.1","ERS537518",
               "GCF-001068415.1","GCF-001069835.1","GCF-000229005.1","GCF-011765545.1","GCF-011765485.1",
               "GCF-011765505.1","GCF-011765525.1","GCF-000522945.1","AHSK01-M2A","GCF-000522565.1",
               "AOBZ01","GCF-000522765.1","GCF-014621675.1","GCF-000229225.1","AHTE01-24",
               "GCF-012642085.1","GCF-000522785.1","GCF-000522645.1","GCF-020529645.1","GCF-020529905.1",
               "GCF-020530625.1","GCF-020530245.1","GCF-901875565.1","GCF-000229085.1","GCF-000229265.1",
               "AHSR01-SM6","GCF-000228985.1","GCF-000229105.1","GCF-000229785.1","AHSD01-M21",
               "GCF-020529305.1","GCF-012642145.1","GCF-020529945.1","GCF-900636835.1","GCF-012641805.1",
               "GCF-012641885.1","GCF-020529365.1","GCF-000229365.1","GCF-000229865.1","AHSJ01-NLML9",
               "GCF-020530305.1","SRS3070067","GCF-000522805.2","GCF-020531255.1","GCF-020530865.1",
               "GCF-020531385.1","GCF-020530205.1","GCF-020531145.1","GCF-020530685.1","GCF-020531605.1",
               "GCF-020529545.1","GCF-020530105.1","GCF-020529785.1","GCF-000091645.1","GCF-020531635.1",
               "GCF-020530405.1","GCF-000229845.1","AHSO01-W6","AHSY01-14D","GCF-002157665.1",
               "GCF-020529495.1","GCF-020530465.1","GCF-020530905.1","AHSC01-G123","GCF-000229525.1",
               "GCF-000229905.1","AHSH01-NLML4","GCF-000522665.1","AHTB01-B","GCF-020530355.1",
               "GCF-006386535.1","GCF-012641625.1","AHRS01-2VS1","GCF-000229565.1","GCF-000230045.1",
               "GCF-000228905.1","AOCA01","GCF-000228745.1","GCF-000229745.1","GCF-000229885.1",
               "GCF-000229545.1","GCF-000228765.1","GCF-012641975.1","GCF-000228925.1","GCF-000228845.1",
               "GCF-012642015.1","GCF-012642165.1","GCF-000229465.1","GCF-000230165.1","GCF-020531345.1",
               "GCF-000229385.1","GCF-000229425.1","GCF-012641565.1","AHSW01-NLML1","AHRQ01-4VF1",
               "GCF-020529275.1","SRS3070105","SRS833653","SRS833655","GCF-020529465.1",
               "GCF-020530285.1","GCF-020531005.1","GCF-020531065.1","GCF-020530345.1","GCF-020530365.1",
               "ERS537441","GCF-020530745.1","GCF-021013185.1","GCF-020529345.1","GCF-020529485.1",
               "GCF-020530445.1","GCF-020531505.1","GCF-020530565.1","GCF-000284575.1","ERS537492",
               "ERS537553","AOCC01","AHSU01-U2A","AHRI01-OMZ175","GCF-002212925.1","GCF-020529825.1",
               "GCF-020530765.1","AHSV01-NLML8","GCF-020530225.1","ERS537485","AHSN01-NV1996",
               "GCF-002212845.1","GCF-008831365.1","ERS537496","GCF-020529385.1","ERS537548","ERS537549",
               "GCF-000522705.1","GCF-020531225.1","SRS3070079","GCF-012642185.1","GCF-020530075.1",
               "GCF-020531465.1","AHSI01-NLML5","GCF-012641845.1","GCF-014842815.1","GCF-014842815.3",
               "GCF-020530805.1","GCF-020531665.1","GCF-020530985.1","GCF-020529625.1","GCF-020531525.1",
               "GCF-002213005.1","GCF-020530165.1","GCF-012642405.1","GCF-020530725.1","GCF-000228965.1",
               "SRS833646","SRS833647","GCF-008831345.1","AHSF01-N34","GCF-018619415.1","AHTD01-SM1",
               "GCF-020530645.1","AHSX01-1ID3","AHTA01-66-2A","GCF-020530125.1","GCF-002212855.1",
               "GCF-002212995.1","GCF-020529565.1","GCF-900459345.1","GCF-012641085.1","GCF-000229345.1",
               "GCF-000229685.1","AHSQ01-SF14","GCF-020530185.1","GCF-020530925.1","GCF-020529685.1",
               "GCF-020529965.1","GCF-020531425.1","GCF-020529525.1","GCF-020530945.1","GCF-020531105.1",
               "AHRR01-15VF2","AHST01-ST6","GCF-000230205.1","AHRC01-S1B","AHRM01-3SN1","GCF-000522605.1",
               "GCF-000229145.1","GCF-000229025.1","GCF-000229765.1","AHRH01-M230","GCF-000522685.1",
               "GCF-012641615.1","GCF-012642225.1","GCF-012641925.1","GCF-020529705.1","AHSP01-SF1",
               "AHTC01-SM4","AHRT01-11VS1","GCF-000229305.1","GCF-002995555.1","GCF-020530665.1",
               "GCF-020531285.1","ERS537533","GCF-020530605.1","ERS537646","ERS537450","ERS537491",
               "GCF-020530045.1","GCF-020531245.1","ERS537476","GCF-008831325.1","GCF-020529765.1",
               "ERS537475","ERS537636","AHRW01-NVAB","AOCB01","GCF-002212905.1","GCF-900475095.1",
               "GCF-000375505.1","GCF-019048645.1","GCF-006739205.1","GCF-012273155.1","AHRG01-R221",
               "GCF-000229185.1","GCF-000522865.1","GCF-012641965.1","AHRN01-2ST1","AHSZ01-21",
               "GCF-018588765.1","GCF-000229725.1","GCF-000230065.1","GCF-012642365.1","GCF-000271865.1",
               "GCF-015670285.1","AHRX01-A9","GCF-020529245.1","GCF-020529405.1","GCF-020530705.1",
               "GCF-020530875.1","AGWE01-U2B","AHRB01-8ID3","AOBY01","GCF-000817065.1",
               "GCF-000229985.1","GCF-000230025.1","AHRV01-NFSM2","GCF-012642035.1","GCF-000230125.1",
               "GCF-000230185.1","GCF-000228825.1","GCF-012642275.1","GCF-020529665.1","AOBX01",
               "GCF-020529445.1","AHSM01-N66","GCF-000229065.1","GCF-000230145.1","GCF-000229165.1",
               "GCF-000229645.1","GCF-000229665.1","AHSA01-A19","AHSB01-U138","AHSS01-ST1",
               "GCF-015668935.1","GCF-020530525.1","GCF-020529725.1","GCF-020531485.1","GCF-012641825.1",
               "GCF-000229125.1","GCF-012641905.1","GCF-012642105.1","GCF-012642125.1","GCF-001625005.1",
               "GCF-012642325.1","GCF-015669655.1","GCF-020531125.1","GCF-000229325.1","GCF-000229625.1",
               "GCF-001703615.1","GCF-020529415.1","GCF-020529795.1","GCF-020530325.1","GCF-020531445.1",
               "GCF-020530485.1","GCF-020531355.1","GCF-000522625.1","GCF-020529585.1","GCF-000522725.1",
               "AHRF01-SF12","GCF-020530425.1","GCF-020531045.1","AHRO01-11A1","GCF-000229585.1",
               "GCF-012641665.1","GCF-000229605.1","GCF-012642055.1","GCF-012641935.1","GCF-012642265.1",
               "GCF-015670115.1","GCF-020530495.1","GCF-020529325.1","GCF-020530065.1","GCF-000229805.1",
               "GCF-000229825.1","GCF-000229705.1","GCF-000228945.1","GCF-000229965.1","GCF-000230085.1",
               "GCF-000230105.1","GCF-000229925.1","GCF-000228865.1","GCF-000229245.1","GCF-000229045.1",
               "GCF-012641775.1","GCF-000522905.1","GCF-000522585.1","GCF-020531585.1","GCF-000522825.1",
               "GCF-001558215.1","GCF-009738105.1","GCF-012641855.1","GCF-012642345.1","GCF-001073145.1",
               "GCF-012642355.1","GCF-002212935.1","AHSL01-N3209","GCF-020529915.1","GCF-012642205.1",
               "SRS2403917","GCF-000522925.1","GCF-012642305.1","GCF-020530145.1","AHRJ01-15JP3",
               "GCF-018588785.1","GCF-000230005.1","GCF-000228885.1","GCF-000229505.1","AHSE01-T4",
               "AHRE01-SA41","GCF-000229205.1","GCF-002212885.1","GCF-000496555.1","GCF-002213035.1",
               "GCF-020530785.1","GCF-020529225.1","GCF-020530585.1","GCF-020529595.1","GCF-020531085.1",
               "GCF-020529865.1","GCF-020531025.1","GCF-020531325.1","GCF-020531205.1","GCF-020531305.1",
               "GCF-000522845.1","GCF-020529255.1","GCF-902365065.1","GCF-020529985.1","GCF-020530845.1",
               "GCF-020530025.1","GCF-002212965.1","GCF-000496535.1","GCF-002213065.1","ERS537424",
               "ERS537571","ERS537484","GCF-002155285.1","GCF-020531405.1")
M <- M[tip_order,]

# get data for annotations
b <- fread(beast_clade_df, data.table=F)
rownames(b) <- b$ID
b$ID <- NULL
anno <- b %>% dplyr::select(Clade)
l <- list(Clade=c(Killuragh = "black", 
                  clade1="#7F3C8D", clade2="#11A579",	clade3="#3969AC",	clade4="#F2B701",	clade5="#E73F74",
                  clade6="#BE0032", clade7="#80BA5A",clade8="#A5AA99"))
library(pheatmap)
pheatmap(M, color=c("white","black"), annotation_row=anno, annotation_colors=l,show_colnames=T, show_rownames=F, cluster_rows=F, cluster_cols=F)
rowSums(M)

#######

## Analysis with tree data
rm(list=ls())


b <- fread(beast_dat, data.table=F)
clades <- fread(beast_clade_df, data.table =F)
mutacins <- fread(mutacin_data, data.table = F)
colnames(mutacins) <- c("label","mutacinI","mutacinII","mutacinIII","mutacinSmb","mutacinK8","mutacinIV","mutacinV")
colnames(clades) <- c("label","Clade")

#merge clade and mutacin info (speed)

merge <- merge(mutacins, clades)

# get ancestral data for each individual

anc <- merge(merge,b)
## Part 1: mutacins & TMRCA
# rename columns and re-merge with b to get ancestor height data

anc <- anc %>% select(c("label","mutacinI","mutacinII","mutacinIII","mutacinSmb","mutacinK8","mutacinIV","mutacinV","Clade","parent"))
colnames(anc) <- c("label","mutacinI","mutacinII","mutacinIII","mutacinSmb","mutacinK8","mutacinIV","mutacinV","Clade","node")

parent <- merge(anc, b, by.x="node", by.y="node")
parent <- parent %>% select("node","label.x","mutacinI","mutacinII","mutacinIII","mutacinSmb","mutacinK8","mutacinIV","mutacinV","Clade", "branch.length",    
                            "height",          
                            "height_0.95_HPD","height_median","height_range","length","length_0.95_HPD","length_median",    
                            "length_range","posterior","rate","rate_0.95_HPD","rate_median","rate_range")


ggplot(parent, aes(x=height_median))  + geom_density() + facet_grid(~Clade) + theme_bw()

ggplot(parent, aes(x=height_median, col=Clade, fill=Clade))  + geom_histogram() + facet_grid(~Clade) + theme_bw() +
  scale_color_manual(values=c( "#7F3C8D", "#11A579",	"#3969AC",	"#F2B701",	"#E73F74","#BE0032", "#80BA5A","#A5AA99","black")) +
  scale_fill_manual(values=c("#7F3C8D", "#11A579",	"#3969AC",	"#F2B701",	"#E73F74","#BE0032", "#80BA5A","#A5AA99","black")) + theme(legend.position = "none")

# Mutacin plots

mrca_mut1 <- parent %>% filter(mutacinI == 1)
mrca_mut2 <- parent %>% filter(mutacinII == 1)
mrca_mut3 <- parent %>% filter(mutacinIII == 1)
mrca_mut4 <- parent %>% filter(mutacinIV == 1)
mrca_mut5 <- parent %>% filter(mutacinV == 1)
mrca_mutsmb <- parent %>% filter(mutacinSmb == 1)
mrca_mutk8 <- parent %>% filter(mutacinK8 == 1)

mrca_mut_neg <- parent %>% filter(mutacinI == 0)%>% filter(mutacinII == 0)%>% filter(mutacinIII == 0)%>% 
  filter(mutacinIV == 0)%>% filter(mutacinV == 0)%>% filter(mutacinSmb == 0)%>% filter(mutacinK8 == 0)

ggplot(parent, aes(x=height_median)) + geom_density(fill="grey", alpha=0.4) #+

ggplot(parent, aes(x=height)) + geom_density(fill="grey", alpha=0.4) #+

ggplot() + geom_density(data=mrca_mut1,aes( x=height_median), fill="#7F3C8D", alpha=0.4) +
  geom_density(data=mrca_mut2,aes( x=height_median), fill="#F97B72", alpha=0.4) +
  geom_density(data=mrca_mut3,aes( x=height_median), fill="#A5AA99", alpha=0.4) +
  geom_density(data=mrca_mut4,aes( x=height_median), fill="#11A579", alpha=0.4) +
  geom_density(data=mrca_mut5,aes( x=height_median), fill="#3969AC", alpha=0.4) +
  geom_density(data=mrca_mutsmb,aes( x=height_median), fill="#F2B701", alpha=0.4) +
  geom_density(data=mrca_mutk8,aes( x=height_median), fill=	"#E73F74", alpha=0.4) +
  theme_bw()

## Part 2 (fig S10) Mutacins vs branch length
# Look at branch length distribution vs mutacin profile
df <-  merge(merge,b)
branch_data <- df %>% select(c("label","mutacinI","mutacinII","mutacinIII","mutacinSmb","mutacinK8","mutacinIV","mutacinV","Clade","node","branch.length"))

# add a column for mutacin positive/negative 
branch_data <- branch_data %>% mutate(mut_status = ifelse((mutacinI == 0 & mutacinII == 0 & mutacinIII == 0 & mutacinIV == 0 & mutacinV == 0 & mutacinSmb == 0 & mutacinK8 == 0),"negative","positive"))
branch_data$mut_status <- factor(branch_data$mut_status, levels = c("positive","negative"))

ggplot(branch_data, aes(x=mut_status, y=branch.length)) + geom_violin() + theme_bw()

library(ggsignif)


ggplot(branch_data, aes(x=mut_status, y=branch.length)) + geom_violin(fill="#ADD8E6") + geom_point() + theme_bw() + 
  geom_signif(comparisons = list(c("positive","negative")), map_signif_level = FALSE, test = "wilcox.test")

ggplot(branch_data, aes(x=mut_status, y=branch.length)) + geom_boxplot()  + theme_bw() + 
  geom_signif(comparisons = list(c("positive","negative")), map_signif_level = FALSE, test = "wilcox.test")


##########################

# Geog structure analysis
rm(list=ls())
dev.off()
# Need to plot tree & associated metadata
tr <- read.nexus(beast_tree)
b <- fread(beast_dat, data.table=F)
clades <- fread(beast_clades_df, data.table =F)
mutacins <- fread(mutacin_data, data.table = F)
colnames(mutacins) <- c("label","mutacinI","mutacinII","mutacinIII","mutacinSmb","mutacinK8","mutacinIV","mutacinV")
colnames(clades) <- c("label","Clade")

#merge clade and mutacin info (speed)

merge <- merge(mutacins, clades)
meta <- fread(metadata_tab, data.table = F)
colnames(meta) <- c("label","Strain","Isolation_Source","Host","Date","Country","Continent")

anc <- merge(meta, b)
 
tr_order <- Preorder(tr)


t1 <- Subtree(tr_order, 812)
t2 <- Subtree(tr_order, 791)
t3 <- Subtree(tr_order, 798)
t4 <- Subtree(tr_order, 688)
t5 <- Subtree(tr_order, 679)
t6 <- Subtree(tr_order, 601)
t7 <- Subtree(tr_order, 415)
t8 <- Subtree(tr_order, 598)

tr_grouped <- groupClade(tr_order, c(812,791,798,688,679,601,415,598))


tree_coloured <- ggtree(tr_grouped, aes(color=group)) + theme(legend.position='none') + coord_cartesian(clip = 'off') + 
  theme_tree2(plot.margin=margin(6, 120, 6, 6)) +
  theme(legend.position=c(.1, .8))  +
  scale_color_manual(values=c("black", "#7F3C8D", "#11A579",	"#3969AC",	"#F2B701",	
                              "#E73F74","#BE0032", "#80BA5A","#A5AA99","#008695",		"#F97B72"))

# put killuragh on top
tree_coloured <- tree_coloured  %>% rotate(412)%>% rotate(686) %>% rotate(687)
revts_plt <- revts(tree_coloured)

###

library(TreeTools)

tip_labels <- tr_order["tip.label"]
anc_list <- ListAncestors(tr_order[["edge"]][, 1],tr_order[["edge"]][, 2], node=NULL)
# MRCA(x1,x2,anc_list) where x1 and x2 are every possible combination of numbers in len(tip_labels)
combos <- expand.grid(seq(1:length(tip_labels$tip.label)),seq(1:length(tip_labels$tip.label)))  %>% filter(Var1 != Var2)
MRCA_df <- combos  %>% rowwise() %>% mutate(mrca=MRCA(Var1, Var2, anc_list))
colnames(MRCA_df) <- c("Var1","Var2","node")
MRCA_df <- MRCA_df %>% mutate(ID1 = tip_labels$tip.label[Var1]) %>% mutate(ID2 = tip_labels$tip.label[Var2])
MRCA_df <- merge(MRCA_df, b, by="node")
MRCA_df <- MRCA_df %>% select("node","ID1","ID2","height","height_median")
# merge so you have branch length data as well
colnames(MRCA_df) <- c("mrca_node","label","ID2","height","height_median")
MRCA_df <- merge(MRCA_df, b, by="label") %>% select(c("mrca_node","label","ID2","height.x","height_median.x", "length", "length_median"))
colnames(MRCA_df) <- c("mrca_node","ID1","label","height","height_median", "length_ID1", "length_median_ID1")
MRCA_df <- merge(MRCA_df, b, by="label") %>% select(c("mrca_node","label","ID1","height.x","height_median.x", "length_ID1", "length_median_ID1","length", "length_median"))
colnames(MRCA_df) <- c("mrca_node","ID1","ID2","height","height_median", "length_ID1", "length_median_ID1","length_ID2", "length_median_ID2")

# merge with country info
country <- meta %>% select(label, Country, Continent)
colnames(country) <- c("ID1", "Country1", "Continent1")
MRCA_df <- merge(MRCA_df, country)
colnames(country) <- c("ID2", "Country2", "Continent2")
MRCA_df <- merge(MRCA_df, country)
#fwrite(MRCA_df, "mrca_df.with_mrca_height_and_branch_lengths.and_country.csv")


head(MRCA_df)
latlong <- fread(lat_long_df)
colnames(latlong) <- c("Country1", "lat1","long1")
MRCA_df <- merge(MRCA_df, latlong)
colnames(latlong) <- c("Country2", "lat2","long2")
MRCA_df <- merge(MRCA_df, latlong)

head(MRCA_df)
library(geosphere)
MRCA_df <- MRCA_df %>% rowwise() %>% mutate(distance = distm(c(long1,lat1),c(long2,lat2)))
#fwrite(MRCA_df, "mrca_df.with_distances.csv")

ggplot(MRCA_df, aes(x=height_median, y=distance)) + geom_point() + geom_smooth(method='lm')
past_100_yrs <- MRCA_df %>% filter(height_median < 100 & height_median > 000)
ggplot(past_100_yrs, aes(x=height_median, y=distance)) + geom_point() + geom_smooth(method='lm')

diff_country <- MRCA_df %>% filter(Country1 != Country2)
ggplot(diff_country, aes(x=height_median, y=distance)) + geom_point() + geom_smooth(method='lm')
past_100_diff <-  diff_country %>% filter(height_median < 100 & height_median > 0)
ggplot(past_100_yrs, aes(x=height_median, y=distance)) + geom_point() + geom_smooth(method='lm')

date_mat <- meta %>% select("label","Country")

cor.test(past_100_yrs$height_median, past_100_yrs$distance, na.action="na.omit")
# similarly for other time bins.
cor.test(MRCA_df$height_median, MRCA_df$distance, na.action="na.omit")
rownames(date_mat) <- date_mat$label

date_mat$label <- NULL


# Plot with all countries
gheatmap(revts_plt, date_mat, offset=8, colnames=FALSE) + scale_x_ggtree()
#  4kya
past_4k <- anc %>% filter(length_median <=4000) %>%  select("label","Country.x")
rownames(past_4k) <- past_4k$label
past_4k$label <- NULL
gheatmap(revts_plt, past_4k, offset=8, colnames=FALSE) + scale_x_ggtree()

# 1.5kya
past_1.5k <- anc %>% filter(length_median <=1500) %>%  select("label","Country.x")
rownames(past_1.5k) <- past_1.5k$label
past_1.5k$label <- NULL
gheatmap(revts_plt, past_1.5k, offset=8, colnames=FALSE) + scale_x_ggtree()

# 1kya
past_1k <- anc %>% filter(length_median <=1000) %>%  select("label","Country.x")
rownames(past_1k) <- past_1k$label
past_1k$label <- NULL
gheatmap(revts_plt, past_1k, offset=8, colnames=FALSE) + scale_x_ggtree()

# 500 ya
past_500 <- anc %>% filter(length_median <=500) %>%  select("label","Country.x")
rownames(past_500) <- past_500$label
past_500$label <- NULL
gheatmap(revts_plt, past_500, offset=8, colnames=FALSE) + scale_x_ggtree()
# 250 years
past_250 <- anc %>% filter(length_median <=250) %>%  select("label","Country.x")
rownames(past_250) <- past_250$label
past_250$label <- NULL
gheatmap(revts_plt, past_250, offset=8, colnames=FALSE) + scale_x_ggtree()

past_100 <- anc %>% filter(length_median <=100) %>%  select("label","Country.x")
rownames(past_100) <- past_100$label
past_100$label <- NULL
gheatmap(revts_plt, past_100, offset=8, colnames=FALSE) + scale_x_ggtree()


