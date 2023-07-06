# Sourcetracker proportions
library(tidyverse)
library(reshape2)
library(cols4all)
sourcetracker_dir <- "./data/results/metagenomics/sourcetracker_species_rarefy_2000_17-06-22/"
setwd(sourcetracker_dir)
st2 <- read.csv('mixing_proportions.txt', sep='\t')

order <- c("KGH1-A","KGH1-E","KGH2-B","KGH2-F","KGH6-A","ANN2-B","PN107","PN112","PN113","SRA62","KK1","CE003","CE004","KD026",
  "KD044","KD045","KD046","KD055","KD057","KD059","HI1","HI2","HS1","HS2","HS3","L","O1","O4","1H04","1H06","1H07","1H13",
  "1H14","2H06","2H07","2H10","2H11","2H17","MDV248","LM_213_T","LM_306_T","LM_308_T","LM_309_T","LM_403_T","LM_406_T",
  "LM_913_T","La_Brana","Ajv70","Gok2","Gok4","Gok5","Gok7","Ire8","Gnie1","Gnie2","Kow10","Kow11","Kow12","Kow13","Kow14",
  "Kow15","Kow16","Kow17","Kow18","Kow19","Kow20","Kow21","Kow22","Kow23","Kow24","Kow25","Kow26","Kow27","Kow28","Kow29",
  "Kow30","Kow31","Kow32","Kow33","Kow34","Kow35","Kow36","Kow37","Kow38","Kow39","Kow40","Kow41","Kow42","Kow43","Kow44",
  "Kow45","Kow46","Kow47","Kow48","Kow49","Kow50","Kow51","Kow52","Kow53","Kow54","Kow55","Kow56","Kow57","Kow58","Kow59",
  "Kow60","Kow61","Kow8","Kow9","Leg1","Leg2","Leg3","Leg4","Leg5","Leg6","Leg7","Leg8","Leg9","Mar1","Mar2","Mar3","Mar4",
  "Mar5","Mar6","Mar7","Mar8","Mas1","Mas10","Mas11","Mas12","Mas13","Mas14","Mas15","Mas16","Mas17","Mas18","Mas19","Mas2",
  "Mas20","Mas21","Mas22","Mas23","Mas24","Mas25","Mas26","Mas27","Mas3","Mas4","Mas5","Mas6","Mas7","Mas8","Mas9","Niem1",
  "Niem10","Niem11","Niem12","Niem13","Niem14","Niem15","Niem16","Niem17","Niem18","Niem19","Niem2","Niem20","Niem21","Niem22",
  "Niem23","Niem24","Niem25","Niem26","Niem27","Niem28","Niem29","Niem3","Niem30","Niem31","Niem32","Niem33","Niem34","Niem35",
  "Niem4","Niem5","Niem6","Niem7","Niem8","Niem9","Sow1","Sow10","Sow11","Sow12","Sow13","Sow14","Sow15","Sow16","Sow17","Sow18",
  "Sow19","Sow2","Sow20","Sow21","Sow3","Sow4","Sow5","Sow6","Sow7","Sow8","Sow9","Abusir1433t","Abusir1513t","Abusir1515t",
  "Abusir1518t","Abusir1521t","Abusir1526t","Abusir1529t","Abusir1534t","Abusir1537t","Abusir1538t","Abusir1545t",
  "Abusir1546t","Abusir1547t","Abusir1550t","Abusir1560t","Abusir1561t","Abusir1563t","Abusir1564t","Abusir1568t","Abusir1569t",
  "Abusir1570t","Abusir1579t","Abusir1580t","Abusir1583t","Abusir1584t","Abusir1587t","Abusir1590t","Abusir1591t","Abusir1595t",
  "Abusir1596t","Abusir1599t","Abusir1601t","Abusir1602t","Abusir1607t","Abusir1612t","Abusir1615t","Abusir1616t","Abusir1617t",
  "Abusir1618t","Abusir1620t","Abusir1621t","Abusir1623t","Abusir1627t","Abusir1631t","Abusir1633t","Abusir1635t","Abusir1636t",
  "Abusir1647t","Abusir1650t","Abusir1651t","Abusir1654t","Abusir1655t","Abusir1660t","Abusir1665t","Abusir1671t","Abusir1672t",
  "Abusir1673t","Abusir3533t","Abusir3536t","Abusir3544t","Abusir3552t","Abusir3578t","Abusir3610t","Abusir1511c","Abusir1519c",
  "Abusir1574c","Abusir1594c","SA-001","SA-002","SA-003","SA-004","SA-005","SA-006","SA-007","SA-008","SA-009","SA-010","SA-011",
  "DRT001","RUV001","SMD017","SMD046","SMD051","WIG001","KT05","KT13","KT26","KT28","KT32","CS01","CS02","CS03","CS04","CS05",
  "CS06","CS07","CS08","CS09","CS10","CS11","CS12","CS13","CS14","CS15","CS16","CS17","CS18","CS19","CS20","CS21","CS22","CS23",
  "CS24","CS26","CS27","CS28","CS29","CS30","CS31","CS32","CS33","CS34","CS35","CS36","CS37","CS38","CS39","CS40","CS42","CS43",
  "CS44","CS45","CS46","CS47","CS48","B61","G12")

meltst <- melt(st2)

meltst$SampleID <- factor(meltst$SampleID, levels=order)
source_colours <- c("#C9DBAA","#7B848F","#FFFF99","#C4A484","#D04D60","#EF6F6A","lightgrey")

colnames(meltst) <- c("Sample_ID","Source","Proportion")

# Sourcetracker plot
# may need to change dimensions
pdf("sourcetracker_plot.pdf")
ggplot(meltst, aes(x=Sample_ID,y=Proportion,fill=Source)) + geom_bar(position="fill",stat="identity") +
  theme_classic() + scale_fill_manual(values=source_colours) +
  ggtitle("Killuragh Sourcetracker") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=4)) 
dev.off()


