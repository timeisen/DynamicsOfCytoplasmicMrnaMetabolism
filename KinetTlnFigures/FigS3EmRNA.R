#FigS3EmRNA fig
#halflife_figure
library(ggplot2)
library(plyr)
library(reshape2)
options("scipen"=5)
source("/lab/solexa_bartel/teisen/RNAseq/Scripts/general/ggplot_theme.R")

actD_samples<-read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/processed_files/halflives/miR-1_actD_halflives_logspace_fits_indiv_baseline.txt",head=TRUE)
# hl_5EU<-read.table("/lab/solexa_bartel/eichhorn/5EU_RPF/miR-1_miR-155_analysis/Halflife_analysis/Halflife_measurements/miR-155_RNA_minusdox_halflife_SS.txt",head=TRUE,sep="\t")
hl_5EU<-read.table("processed_files/halflifeComparisons/miR-155_minu_PAL_halflives_logspace_global_offset_SS_RNAseqScaling.txt",head=TRUE,sep="\t")

# hl_5EU$beta <- log(2)/hl_5EU$halflife

hl<-merge(actD_samples,hl_5EU,by="accession")

# hl<-hl[which(hl$residual_actD<10E-11),]

hl[hl$halflife_t < 0.1,]$halflife_t = 0.1
hl[hl$halflife_t > 1E2,]$halflife_t = 1E2
hl[hl$halflife_actD < 0.1,]$halflife_actD = 0.1
hl[hl$halflife_actD > 1E2,]$halflife_actD = 1E2

p1<-ggplot(hl,aes(x=halflife_t,y=halflife_actD)) +
	geom_point(alpha=0.2,size=0.5,shape = 16) +
	theme_tim() +
	scale_x_continuous(name=NULL,trans="log10",limits=c(0.1,1E2),labels=log_ticks(0.1,1E2)[[1]],breaks=log_ticks(0.1,1E2)[[2]],expand=c(0,0)) +
	scale_y_continuous(name=NULL,trans="log10",limits=c(0.1,1E2),labels=log_ticks(0.1,1E2)[[1]],breaks=log_ticks(0.1,1E2)[[2]],expand=c(0,0)) +
	geom_abline(linetype="dashed",color="grey",slope=1,intercept=0)

print(nrow(hl))
print(cor(hl$halflife_actD,log(2)/hl$beta_t/60,method="spearman"))
print(median(hl$halflife_t))
print(median(hl$halflife_actD))


ggsave(p1,file="/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/figures/version_9/FigS6C.pdf",width=2,height=2)