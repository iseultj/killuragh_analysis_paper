rm(list=ls())
wd <- './data/results/denticola/heteroplasmy/'

setwd(wd)
library(tidyverse)
library(data.table)

sums_noncompaln <- read.csv('aggregated_maj_minor_sums.non_comp_aln.ancients_only.csv')
sums_ca <- read.csv('aggregated_maj_minor_sums.comp_aln.22-3-23.csv')
sums_old_ca <- read.csv('aggregated_maj_minor_sums.refseq_comp_aln.csv')
sums_noncompaln$X <- NULL
sums_old_ca$X <- NULL
sums_ca$X <- NULL
# sort Position

sums_noncompaln <- sums_noncompaln[order(sums_noncompaln$Position),]
sums_old_ca <- sums_old_ca[order(sums_old_ca$Position),]
sums_ca <- sums_ca[order(sums_ca$Position),]

# sliding windows


rolling_summary <- function(DF, time_col, fun, window_size, step_size, min_window=min(DF$time_col)) {
  # time_col is name of time (position) column
  # fun is function to apply to the subsetted data frames
  # min_window is the start time of the earliest window
  
  times <- DF[, time_col]
  m_t <- max(times) 
  # window_starts is a vector of the windows' minimum times
  window_starts <- seq(from=min_window, to=m_t, by=step_size)
  
  # The i-th element of window_rows is a vector that tells us the row numbers of
  # the data-frame rows that are present in window i 
  window_rows <- lapply(window_starts, function(x) { which(times>=x & times<x+window_size) })
  
  window_summaries <- sapply(window_rows, function(w_r) fun(DF[w_r, ]))
  #print(window_summaries)
  data.frame(start_pos=window_starts, end_pos=window_starts+window_size, summary=window_summaries)
}


het <- function(df) {
  major_sums <- sum(df$major_sums)
  minor_sums <- sum(df$minor_sums)
  heteroplasmy <- minor_sums/major_sums
  return(heteroplasmy)
}

het_no_md <- function(df) {
  major_sums <- sum(df$major_sums)
  minor_sums <- sum(df$minor_sums_no_md)
  heteroplasmy <- minor_sums/major_sums
  return(heteroplasmy)
}

rolling_10kb_window_1kb_step <- rolling_summary(sums_old_ca,"Position",fun=het,window_size=10000, step_size=1000, min_window=220)
rolling_10kb_window_1kb_step_no_compaln <- rolling_summary(sums_noncompaln,"Position",fun=het,window_size=10000, step_size=1000, min_window=220)
rolling_10kb_window_1kb_step_new_compaln <- rolling_summary(sums_ca,"Position",fun=het,window_size=10000, step_size=1000, min_window=220)
colnames(rolling_10kb_window_1kb_step) <- c("start","end","het")
colnames(rolling_10kb_window_1kb_step_new_compaln) <- c("start","end","het")
colnames(rolling_10kb_window_1kb_step_no_compaln) <- c("start","end","het")
rolling_10kb_window_1kb_step_no_md <- rolling_summary(sums_old_ca,"Position",fun=het_no_md,window_size=10000, step_size=1000, min_window=220)
rolling_10kb_window_1kb_step_no_compaln_no_md <- rolling_summary(sums_noncompaln,"Position",fun=het_no_md,window_size=10000, step_size=1000, min_window=220)
rolling_10kb_window_1kb_step_new_compaln_no_md <- rolling_summary(sums_ca,"Position",fun=het_no_md,window_size=10000, step_size=1000, min_window=220)
colnames(rolling_10kb_window_1kb_step_no_md) <- c("start","end","het")
colnames(rolling_10kb_window_1kb_step_no_compaln_no_md) <- c("start","end","het")
colnames(rolling_10kb_window_1kb_step_new_compaln_no_md)<- c("start","end","het")
c1 <- ggplot(rolling_10kb_window_1kb_step, aes(x=start, y=het)) + geom_point(fill="#537d90", col="#537d90", alpha=0.7) + theme_bw() + ggtitle("10kb windows, 1kb step, competitive alignment") + ylim(0,0.25)
c2 <- ggplot(rolling_10kb_window_1kb_step_no_compaln, aes(x=start, y=het)) + geom_point(fill="#7F3C8D", col="#7F3C8D", alpha=0.7) + theme_bw() + ggtitle("10kb windows, 1kb step, normal alignment")  + ylim(0,0.25)
c4 <- ggplot(rolling_10kb_window_1kb_step_no_md, aes(x=start, y=het)) + geom_point(fill="#537d90", col="#537d90", alpha=0.7) + theme_bw() + ggtitle("10kb windows, 1kb step, competitive alignment, no MD") + ylim(0,0.25)
c5 <- ggplot(rolling_10kb_window_1kb_step_no_compaln_no_md, aes(x=start, y=het)) + geom_point(fill="#7F3C8D", col="#7F3C8D", alpha=0.7) + theme_bw() + ggtitle("10kb windows, 1kb step, normal alignment, no MD")  + ylim(0,0.25)
c3 <- ggplot(rolling_10kb_window_1kb_step_new_compaln, aes(x=start, y=het)) + geom_point(fill="#BE0032", col="#BE0032", alpha=0.7) + theme_bw() + ggtitle("10kb windows, 1kb step, new competitive alignment") + ylim(0,0.25)
c6 <- ggplot(rolling_10kb_window_1kb_step_new_compaln_no_md, aes(x=start, y=het)) + geom_point(fill="#BE0032", col="#BE0032", alpha=0.7) + theme_bw() + ggtitle("10kb windows, 1kb step, new competitive alignment,no MD") + ylim(0,0.25)
c3
library(cowplot)
plot_grid(c1,c2,c3,c4,c5,c6)
rolling_1kb_window_1kb_step_new_compaln <- rolling_summary(sums_ca,"Position",fun=het,window_size=1000, step_size=1000, min_window=220)
colnames(rolling_1kb_window_1kb_step_new_compaln) <- c("start","end","het")
ggplot(rolling_1kb_window_1kb_step_new_compaln, aes(x=start, y=het)) + geom_point(fill="#BE0032", col="#BE0032", alpha=0.7) + theme_bw() 
fwrite(rolling_1kb_window_1kb_step_new_compaln,"rolling_1kb_window_1kb_step.heteroplasmy.competitive_alignment.csv")
