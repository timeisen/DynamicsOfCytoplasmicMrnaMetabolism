#No longer a figure, was meant originally to investigate individual genes tail length and expression relationship.

library(tidyverse)
library(data.table)
source("/lab/solexa_bartel/teisen/RNAseq/Scripts/general/ggplot_theme.R")
annot<-read.table("/lab/solexa_bartel/teisen/RNAseq/Annotation_files/TAIL_seq/symbol_to_accession.txt",head=TRUE)
rate_constants <- read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/rate_constant_measurements/miR-155_minus_samples_v74_plogis_linked_run1.txt",head=1)
colnames(rate_constants)<-c("accession","st","a","k","r")
#need to include "Atf4"
genes <- c("Actb","Eef2","Cyr61","Plk2","Atf4","Slc3a2","Thbs1","Ctgf")

hl<-read.table("/lab/solexa_bartel/eichhorn/5EU_RPF/miR-1_miR-155_analysis/Halflife_analysis/Halflife_measurements/miR-155_RNA_minusdox_halflife_SS.txt",head=TRUE,sep="\t")
hl <- hl[,c(1,3)]

TL <- read.table("model_input_files/miR-155_minus_sample_mean_tails_50tags_with_expr.txt") ##RAW DATA, not background subtracted
TL <- read.table("model_input_files/miR-155_minus_sample_mean_tails_50tags_with_expr_background_subtracted_v5.txt") ##Background subtracted, not background subtracted

# actb_tl <- c(
# 	222.96,
# 	226.07,
# 	196.4,
# 	193.6,
# 	174.85,
# 	159.86)

# actb_tl <- c(140.73,
# 142.17,
# 128.22,
# 126.88,
# 117.75,
# 110.26)

# actb_expr <- c(
# 	142.68/1356.56,
# 	145.14/1059.97,
# 	63.97/218.07,
# 	86.32/99.95,
# 	99.81/57.97,
# 	135.54/50.94)

# eef2_tl <- c(
# 	250.12,
# 	226.18,
# 	195.73,
# 	169.75,
# 	163.67,
# 	134.58)

# eef2_tl <- c(156.74,
# 145.72,
# 131.22,
# 118.35,
# 115.28,
# 100.1)



# eef2_expr <- c( 
# 	21.15/1356.56,
# 	26.16/1059.97,
# 	20/218.07,
# 	18.85/99.95,
# 	17.94/57.97,
# 	43.13/50.94)

colnames(TL) <- c("accession","t40m","t1h","t2h","t4h","t8h","tSs","e40m","e1h","e2h","e4h","e8h","eSs")

m <- merge(TL,hl)
m <- merge(m,annot)

ms <- m[which(m$symbol %in% genes),]#subset

mst <- merge(ms,annot)
mst <- merge(mst,rate_constants)
print(mst[,c(1,2,13,14,15:19)])

mdf <- data.frame(
	accession = rep(ms$accession,each = 6),
	symbol = rep(ms$symbol,each=6),
	tl = c(t(as.matrix(ms[,2:7]))),
	expr = c(t(as.matrix(ms[,8:13]))),
	tp = rep(c(40,60,120,240,480,6000),nrow(ms))
	)

# mdf <- cbind(mdf,c(actb_tl,eef2_tl),c(actb_expr,eef2_expr))
# colnames(mdf)[c(6,7)]<-c("blot_tl","blot_expr")

p6<-ggplot(mdf,aes(x=tl,y=expr,color=as.factor(tp)))+geom_path(alpha=.5,size=1,aes(group=accession))+facet_wrap(~symbol,nrow=1)+scale_y_continuous(trans="log10",name=NULL,expand=c(0,0),limits=c(1E4,1E10))+theme_tim()+scale_x_continuous(name="Tail Length (nt)",expand=c(0,0),limits=c(0,250))+guides(colour = guide_legend(override.aes = list(alpha=1)))+geom_point(size=0.1)#+scale_color_manual(values=cols)+theme()

p7<-ggplot(mdf,aes(x=tp,y=tl,color=as.factor(tp)))+geom_path(alpha=.5,size=1,aes(group=accession))+facet_wrap(~symbol,ncol=6)+scale_y_continuous(name=NULL,expand=c(0,0))+theme_tim()+scale_x_continuous(name="Tail Length (nt)",expand=c(0,0))+guides(colour = guide_legend(override.aes = list(alpha=1)))+geom_point(size=0.1)#+scale_color_manual(values=cols)+theme()

ggsave(plot=p6,file="/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/figures/other/expr_tl_pulse.pdf",width=10,height=8,useDingbats=FALSE)

# p8 <- ggplot(mdf,aes(x=tl,y=blot_tl))+geom_point()+facet_wrap(~symbol)+geom_abline(slope=1,intercept=0,linetype="dashed",color="grey")
# p9 <- ggplot(mdf,aes(x=expr,y=blot_expr))+geom_point()+facet_wrap(~symbol)

# ggsave(plot=p8,file="/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/figures/other/tl_vs_blot.pdf",width=5,height=5,useDingbats=FALSE)
# ggsave(plot=p9,file="/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/figures/other/expr_vs_blot.pdf",width=5,height=5,useDingbats=FALSE)
