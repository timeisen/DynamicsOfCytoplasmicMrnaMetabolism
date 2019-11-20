#This code is for histogram of the half lives, for fig S3B
#This code is for violin plots of the half lives, for fig S3C

library(tidyverse)
library(plyr)
library(reshape2)
source("/lab/solexa_bartel/teisen/RNAseq/Scripts/general/ggplot_theme.R")
IEGs<-read.table("/lab/solexa_bartel/eichhorn/5EU_RPF/miR-1_miR-155_analysis/Halflife_analysis/Halflife_measurements/miR-155_RNA_minusdox_alpha_hl_genenames_iegs.txt")
RPGs<-read.table("/lab/solexa_bartel/eichhorn/5EU_RPF/miR-1_miR-155_analysis/Halflife_analysis/Halflife_measurements/miR-155_RNA_minusdox_alpha_hl_genenames_rgps.txt")
halflife<-read.table("processed_files/halflifeComparisons/miR-155_minu_PAL_halflives_logspace_global_offset_SS_RNAseqScaling.txt",head=TRUE,sep="\t")

halflife<-halflife[,c(1,5)]
colnames(halflife)<-c("accession","halflife")
print(nrow(halflife))
print(median(halflife$halflife))
IEGs<-IEGs[,c(2,4)]
RPGs<-RPGs[,c(2,4)]
colnames(IEGs)<-c("accession","halflife")
colnames(RPGs)<-c("accession","halflife")

# print(length(halflife[which(halflife$halflife<0.001),]$halflife))
halflife[which(halflife$halflife>100),]$halflife = 100
# print(halflife[which(halflife$halflife>100),]$halflife)

halflife$classification<-"All genes"
IEGs$classification<-"IEGs"
RPGs$classification<-"RPGs"

halflife<-rbind(halflife,RPGs,IEGs)

halflife$classification<-factor(halflife$classification,levels=c("All genes","RPGs","IEGs"))
halflife[halflife$halflife < 0.1,]$halflife = 0.1

p1<-ggplot(halflife,aes(x=classification,y=halflife)) + 
	geom_violin(color="black") + 
	theme_tim() + 
	scale_y_continuous(trans="log10",breaks=log_ticks(.1,100)[[2]],labels=log_ticks(.1,100)[[1]],limits=c(0.1,101),name=NULL) + 
	scale_x_discrete(name=NULL)

ggsave(plot=p1,file="/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/figures/version_9/Fig3B.pdf",width=1.5,height=2,useDingbats=FALSE)
#histogram
p2 <- ggplot(halflife,aes(x=halflife,y=..density../50)) + 
	stat_bin(geom="step",bins=50,position='identity',color="black") + 
	theme_tim() + 
	scale_x_continuous(trans="log10",breaks=log_ticks(.1,100)[[2]],labels=log_ticks(.1,100)[[1]],limits=c(0.09,100),name=NULL) + 
	scale_y_continuous(name=NULL,expand=c(0,0))

ggsave(plot=p2,file="/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/figures/version_9/Fig3A.pdf",width=2,height=2,useDingbats=FALSE)
