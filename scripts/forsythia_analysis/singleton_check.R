rm(list=ls())
wd <- './data/intermediate/forsythia/singleton_analysis'
setwd(wd)

library(tidyverse)
library(data.table)

dat <- fread('aggregated_data_all_samples.with_1x.tsv', data.table=F)
dat$filter_sites <- paste0(dat$snpcall_approach,"_",dat$filter)
meta <- fread('UDG_cov_data.tsv', data.table = F)
dat1 <- dat %>% select(ID, filter_sites,n_singletons)
D <- reshape2::dcast(dat1, ID~filter_sites)

# drop KGH not separated

D <- D %>% filter(ID != "KGH")
# drop dog data
D <- D %>% filter(ID != "OH2617-COT023")
# set NA to 0
D[is.na(D)] <- 0

colnames(D)
colnames(meta)
D <- merge(D, meta)

D$tv_ratio_2x <- D$strict_2x_tvs_singletons/D$relaxed_2x_tvs_singletons
D$tv_tn_ratio_2x <- D$strict_2x_tvs_tns_singletons/D$relaxed_2x_tvs_tns_singletons

# set any infinite values to 100
D[sapply(D, is.infinite)] <- 100
library(ggrepel)
ggplot(D, aes(x=tv_ratio_2x, y=tv_tn_ratio_2x, label=ID, col=Qualimap_2x_pc, fill=Qualimap_2x_pc)) + 
  geom_point() + geom_abline() + theme_minimal() + geom_text_repel()

ggplot(D, aes(x=tv_ratio_2x, y=tv_tn_ratio_2x, label=ID, col=UDG, fill=UDG)) + 
  geom_point() + geom_abline() + theme_minimal() + geom_text_repel()

ggplot(D, aes(x=strict_2x_tvs_singletons, y=strict_2x_tvs_tns_singletons, label=ID, col=UDG, fill=UDG)) + 
  geom_point() + geom_abline() + theme_minimal() + geom_text_repel()
ggplot(D, aes(x=strict_1x_tvs_singletons, y=strict_1x_tvs_tns_singletons, label=ID, col=UDG, fill=UDG)) + 
  geom_point() + geom_abline() + theme_minimal() + geom_text_repel()
ggplot(D, aes(x=relaxed_1x_tvs_singletons, y=relaxed_1x_tvs_tns_singletons, label=ID, col=UDG, fill=UDG)) + 
  geom_point() + geom_abline() + theme_minimal() + geom_text_repel()

ggplot(D, aes(x=relaxed_2x_tvs_singletons, y=relaxed_2x_tvs_tns_singletons, label=ID)) + 
  geom_point(aes(col=Qualimap_2x_pc, fill=Qualimap_2x_pc)) + geom_smooth(method="lm") + theme_minimal() + geom_text_repel()

ggplot(D, aes(x=strict_2x_tvs_singletons, y=strict_2x_tvs_tns_singletons, label=ID)) + 
  geom_point(aes(col=Qualimap_2x_pc, fill=Qualimap_2x_pc)) + geom_smooth(method="lm") + theme_minimal() + geom_text_repel()

ggplot(D, aes(x=strict_3x_tvs_singletons, y=strict_3x_tvs_tns_singletons, label=ID)) + 
  geom_point(aes(col=Qualimap_2x_pc, fill=Qualimap_2x_pc)) + geom_smooth(method="lm") + theme_minimal() + geom_text_repel()
ggplot(D, aes(x=relaxed_3x_tvs_singletons, y=relaxed_3x_tvs_tns_singletons, label=ID)) + 
  geom_point(aes(col=Qualimap_2x_pc, fill=Qualimap_2x_pc)) + geom_smooth(method="lm") + theme_minimal() + geom_text_repel()


ggplot(D, aes(x=relaxed_3x_tvs_tns_singletons, y=relaxed_2x_tvs_tns_singletons, label=ID)) + 
  geom_point(aes(col=Qualimap_2x_pc, fill=Qualimap_2x_pc)) + geom_smooth(method="lm") + theme_minimal() + geom_text_repel()

ggplot(D, aes(x=relaxed_3x_tvs_singletons, y=relaxed_2x_tvs_singletons, label=ID)) + 
  geom_point(aes(col=Qualimap_2x_pc, fill=Qualimap_2x_pc)) + geom_smooth(method="lm") + theme_minimal() + geom_text_repel()

ggplot(D, aes(x=strict_2x_tvs_tns_singletons, y=strict_1x_tvs_tns_singletons, label=ID)) + 
  geom_point(aes(col=Qualimap_2x_pc, fill=Qualimap_2x_pc)) + geom_smooth(method="lm") + theme_minimal() + geom_text_repel()

strict_v_relaxed_2x_all <- ggplot(D, aes(x=strict_2x_tvs_tns_singletons, y=relaxed_2x_tvs_tns_singletons, label=ID)) + 
  geom_point(aes(col=Qualimap_2x_pc, fill=Qualimap_2x_pc))  + theme_minimal() + geom_text_repel() + theme(legend.position="none")

strict_v_relaxed_2x_tv <-ggplot(D, aes(x=strict_2x_tvs_singletons, y=relaxed_2x_tvs_singletons, label=ID)) + 
  geom_point(aes(col=Qualimap_2x_pc, fill=Qualimap_2x_pc))  + theme_minimal() + geom_text_repel() + theme(legend.position="none")

strict_v_relaxed_3x_all <- ggplot(D, aes(x=strict_3x_tvs_tns_singletons, y=relaxed_3x_tvs_tns_singletons, label=ID)) + 
  geom_point(aes(col=Qualimap_2x_pc, fill=Qualimap_2x_pc))  + theme_minimal() + geom_text_repel()+ theme(legend.position="none")

strict_v_relaxed_3x_tv <-ggplot(D, aes(x=strict_3x_tvs_singletons, y=relaxed_3x_tvs_singletons, label=ID)) + 
  geom_point(aes(col=Qualimap_2x_pc, fill=Qualimap_2x_pc))  + theme_minimal() + geom_text_repel()+ theme(legend.position="none")

strict_v_relaxed_1x_tv <-ggplot(D, aes(x=strict_1x_tvs_singletons, y=relaxed_1x_tvs_singletons, label=ID)) + 
  geom_point(aes(col=Qualimap_2x_pc, fill=Qualimap_2x_pc))  + theme_minimal() + geom_text_repel()+ theme(legend.position="none")

strict_v_relaxed_1x_all <- ggplot(D, aes(x=strict_1x_tvs_tns_singletons, y=relaxed_1x_tvs_tns_singletons, label=ID)) + 
  geom_point(aes(col=Qualimap_2x_pc, fill=Qualimap_2x_pc))  + theme_minimal() + geom_text_repel()+ theme(legend.position="none")
library(cowplot)
plot_grid(strict_v_relaxed_1x_all,strict_v_relaxed_1x_tv,strict_v_relaxed_2x_all,
          strict_v_relaxed_2x_tv,strict_v_relaxed_3x_all,strict_v_relaxed_3x_tv)

ggplot(dat1, aes(x=ID,y=n_singletons, col=filter_sites,fill=filter_sites)) + geom_point()

# dat1 filter

dat2 <- dat1 %>%  filter(ID != "KGH") %>% filter(ID != "OH2617-COT023")
dat2[is.na(dat2)] <- 0
# set NA to 0

dat2 <- merge(dat2,meta)
dat2$norm_singletons <- dat2$n_singletons/dat2$Qualimap_2x_pc
ggplot(dat2, aes(x=ID,y=n_singletons, col=filter_sites,fill=filter_sites)) + geom_point() + 
  theme_minimal() +
  scale_colour_manual(values=cols12) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

cols12 <- c("#A6CEE3", "#1F78B4","#B2DF8A","#33A02C" ,"#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#CAB2D6" ,"#6A3D9A","#FFFF99","#B15928")
asc <- dat2 %>% filter(filter_sites %in% c("strict_2x_tvs_singletons","strict_2x_tvs_tns_singletons"))
ggplot(dat2, aes(x=ID,y=norm_singletons, col=filter_sites,fill=filter_sites)) + geom_point(alpha=0.7) + 
  scale_fill_manual(values=cols12) + scale_color_manual(values=cols12) + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(asc, aes(x=ID,y=norm_singletons, col=filter_sites,fill=filter_sites)) + geom_point(alpha=0.7) + 
  scale_fill_manual(values=c("#1F78B4","#FB9A99")) + scale_color_manual(values=c("#1F78B4","#FB9A99")) + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
