#Figure 5D, mRNA paper
#Decapping rate histogram.

library(ggplot2)
options("scipen"=5)
source("/lab/solexa_bartel/teisen/RNAseq/Scripts/general/ggplot_theme.R")
options(warn=1)
# hl <- read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/processed_files/halflifeComparisons/miR-155_minu_PAL_halflives_logspace_global_offset_SS_RNAseqScaling.txt",head=TRUE,sep="\t")
o8<-read.table("rate_constant_measurements/miR-155_minus_samples_UNLINKV3_V81H5_run3_high_precision_reformat.txt")
annot<-read.table("/lab/solexa_bartel/teisen/RNAseq/Annotation_files/TAIL_seq/symbol_to_accession.txt",head=TRUE)
o8<-o8[complete.cases(o8),]



colnames(o8) <- c("accession","st_ul","a_ul","k_ul","b_ul","r_ul")
# m<-merge(hl,o8,by="accession")
m<-merge(o8,annot,by="accession")

print(nrow(m))
pLogEl = plogis(230, loc = 262.96905,scale = 11.07612)

print(m[which(m$k_ul < 0.03),])
m$b_ul <- m$b_ul*pLogEl

m[which(m$k_ul > 30),]$k_ul = 30
m[which(m$k_ul < 0.03),]$k_ul = 0.03

m[which(m$b_ul > 3),]$b_ul = 3
print(nrow(m[which(m$b_ul < 0.003),]))
m[which(m$b_ul < 0.003),]$b_ul = 0.003

p1<-ggplot(m,aes(x = k_ul, y = ..count../nrow(m)))+
	stat_bin(geom="step",bins=50,position='identity') +
	scale_x_continuous(name=NULL,trans="log10",breaks=log_ticks(0.01,1000)[[2]],labels=log_ticks(0.01,1000)[[1]]) +
	coord_cartesian(xlim = c(0.03,30)) +
	scale_y_continuous(expand=c(0,0)) +
	theme_tim() +
	theme(axis.line.x = element_blank()) + geom_segment(aes(x = 0.03, xend = 30, y = 0, yend = 0))

p2<-ggplot(m,aes(x = b_ul, y = ..count../nrow(m)))+
	stat_bin(geom="step",bins=50,position='identity') +
	scale_x_continuous(name=NULL,trans="log10",breaks=log_ticks(0.001,1000)[[2]],labels=log_ticks(0.001,1000)[[1]]) +
	coord_cartesian(xlim = c(0.003,3)) +
	scale_y_continuous(expand=c(0,0)) +
	theme_tim() +
	theme(axis.line.x = element_blank()) + geom_segment(aes(x = 0.003, xend = 3, y = 0, yend = 0))

ggsave(plot=p1,file="/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/figures/version_9/Fig5Ahist.pdf",width=2,height=2)
ggsave(plot=p2,file="/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/figures/version_9/Fig5Bhist.pdf",width=2,height=2)
