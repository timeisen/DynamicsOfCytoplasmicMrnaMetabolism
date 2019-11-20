#Fig 1CD miRNA paper, histograms. 

library(plyr)
library(reshape2)
library(tidyverse)
library(spatstat)
source("ggplot_theme.R")

#Histograms of miRNAs.
# args <- commandArgs(trailingOnly=TRUE)

#This code was modified on 2019 02 27 to create CDFs from direct ligation data. 
#It uses length matched no site data in the following folder: /lab/solexa_bartel/teisen/Tail-seq/PalSeqDirectLig/Analysis/miR_analysis/tpUTRlenMatching

cdf<-function(x){
	cdf = sapply(1:length(x),function(n){sum(x[1:n])})
	return(cdf/cdf[length(cdf)])}

closest_number <- function(num,vec){
	return(which(abs(vec-num)==min(abs(vec-num)))[1]) #note added the 1 in case of ties
}

differences <- function(x,m){
	# if(x[1] == 229){
	# print(x[3])
	# print(x[1])
	# print(m[closest_number(x[3],m$plus_nosite),1])
	# print(m[closest_number(x[3],m$plus_nosite),1] - x[1])
	# break
	# }
	return(c(
	m[closest_number(x[3],m$minu_nosite),1] - x[1],
	m[closest_number(x[3],m$minu_wisite),1] - x[1],
	m[closest_number(x[3],m$plus_nosite),1] - x[1]
	))

	}

args1<-c("/archive/bartel/2015.03.04-11263/solexa_bartel/eichhorn/Ribosome_profiling_paper/Sites/Mouse/miR-1/mmu-miR-1_tpUTR_oneormore.txt","/lab/solexa_bartel/teisen/Tail-seq/PalSeqDirectLig/Analysis/miR_analysis/tpUTRlenMatching/nosite_miR1_tpUTR_20190725miR.txt","directLig_1_site.pdf","1")

args2<-c("/archive/bartel/2015.03.04-11263/solexa_bartel/eichhorn/Ribosome_profiling_paper/Sites/Mouse/miR-155/mmu-miR-155_tpUTR_oneormore.txt","/lab/solexa_bartel/teisen/Tail-seq/PalSeqDirectLig/Analysis/miR_analysis/tpUTRlenMatching/nosite_miR155_tpUTR_20190725miR.txt","directLig_155_site.pdf","155")

args3<-c("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/miR_analysis/miR-1_3T3_25pct_down.txt","/lab/solexa_bartel/teisen/Tail-seq/PalSeqDirectLig/Analysis/miR_analysis/tpUTRlenMatching/nosite_miR1_topsite_20190725miR.txt","directLig_1_topsite.pdf","1")

args4<-c("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/miR_analysis/miR-155_3T3_25pct_down.txt","/lab/solexa_bartel/teisen/Tail-seq/PalSeqDirectLig/Analysis/miR_analysis/tpUTRlenMatching/nosite_miR155_topsite_20190725miR.txt","directLig_155_topsite.pdf","155")

for(args in list(args1,args2,args3,args4)){
# for(args in list(args1)){
	
	if(args[4]=="155"){
	PAL_minu<-read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/single_tag_files/background_subtracted_single_tag_files/SS_minus_norm_miR-155_v7_HYBRID20190725miR.txt",sep="\t",head=FALSE)
	PAL_plus<-read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/single_tag_files/background_subtracted_single_tag_files/SS_plus_norm_miR-155_v7_HYBRID20190725miR.txt",sep="\t",head=FALSE)
	}
	
	if(args[4]=="1"){
	PAL_minu<-read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/single_tag_files/background_subtracted_single_tag_files/SS_minus_norm_miR-1_v7_HYBRID20190725miR.txt",sep="\t",head=FALSE)
	PAL_plus<-read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/single_tag_files/background_subtracted_single_tag_files/SS_plus_norm_miR-1_v7_HYBRID20190725miR.txt",sep="\t",head=FALSE)
	}

	PAL_minu <- PAL_minu[,c(252,1:251)]
	PAL_plus <- PAL_plus[,c(252,1:251)]
	colnames(PAL_minu) <- c("accession",paste0("tl",0:250))
	colnames(PAL_plus) <- c("accession",paste0("tl",0:250))


	wisite<-read.table(args[1])
	nosite<-read.table(args[2])
	
	colnames(nosite)[1]<-"accession"
	colnames(wisite)[1]<-"accession"
	
	colnames(PAL_minu)[1]<-"accession"
	colnames(PAL_plus)[1]<-"accession"
	
	PAL_minu_nosite<-PAL_minu[which(PAL_minu$accession %in% nosite$accession),]
	PAL_plus_wisite<-PAL_plus[which(PAL_plus$accession %in% wisite$accession),]
	PAL_minu_wisite<-PAL_minu[which(PAL_minu$accession %in% wisite$accession),]
	PAL_plus_nosite<-PAL_plus[which(PAL_plus$accession %in% nosite$accession),]
	
	PAL_minu_nosite_median<-median(rowSums(PAL_minu_nosite[-1]))
	PAL_plus_wisite_median<-median(rowSums(PAL_plus_wisite[-1]))
	PAL_minu_wisite_median<-median(rowSums(PAL_minu_wisite[-1]))
	PAL_plus_nosite_median<-median(rowSums(PAL_plus_nosite[-1]))
	
	
	PAL_minu_nosite_high<-which(rowSums(PAL_minu_nosite[-1])>PAL_minu_nosite_median)
	PAL_plus_wisite_high<-which(rowSums(PAL_plus_wisite[-1])>PAL_plus_wisite_median)
	PAL_minu_wisite_high<-which(rowSums(PAL_minu_wisite[-1])>PAL_minu_wisite_median)
	PAL_plus_nosite_high<-which(rowSums(PAL_plus_nosite[-1])>PAL_plus_nosite_median)
	
	#Normalization of the rows that have values above the median value. 
	PAL_minu_nosite[PAL_minu_nosite_high,-1] = PAL_minu_nosite[PAL_minu_nosite_high,-1]/rowSums(PAL_minu_nosite[PAL_minu_nosite_high,-1])*	PAL_minu_nosite_median
	PAL_plus_wisite[PAL_plus_wisite_high,-1] = PAL_plus_wisite[PAL_plus_wisite_high,-1]/rowSums(PAL_plus_wisite[PAL_plus_wisite_high,-1])*	PAL_plus_wisite_median
	PAL_minu_wisite[PAL_minu_wisite_high,-1] = PAL_minu_wisite[PAL_minu_wisite_high,-1]/rowSums(PAL_minu_wisite[PAL_minu_wisite_high,-1])*	PAL_minu_wisite_median
	PAL_plus_nosite[PAL_plus_nosite_high,-1] = PAL_plus_nosite[PAL_plus_nosite_high,-1]/rowSums(PAL_plus_nosite[PAL_plus_nosite_high,-1])*	PAL_plus_nosite_median
	
	m <- data.frame(
		"tail_length" = 0:250,
		"minu_nosite" = colSums(PAL_minu_nosite[,-1])/sum(colSums(PAL_minu_nosite[,-1])),
		"plus_wisite" = colSums(PAL_plus_wisite[,-1])/sum(colSums(PAL_plus_wisite[,-1])),
		"minu_wisite" = colSums(PAL_minu_wisite[,-1])/sum(colSums(PAL_minu_wisite[,-1])),
		"plus_nosite" = colSums(PAL_plus_nosite[,-1])/sum(colSums(PAL_plus_nosite[,-1])))
	
	repSize = 10000
	mSample <- data.frame(
		"tail_length" = 0:250,
		"minu_nosite" = rmultinom(n = 1, size = repSize, prob = colSums(PAL_minu_nosite[,-1])),
		"plus_wisite" = rmultinom(n = 1, size = repSize, prob = colSums(PAL_plus_wisite[,-1])),
		"minu_wisite" = rmultinom(n = 1, size = repSize, prob = colSums(PAL_minu_wisite[,-1])),
		"plus_nosite" = rmultinom(n = 1, size = repSize, prob = colSums(PAL_plus_nosite[,-1])))

	mDataSim <- data.frame(
		"minu_nosite" = rep(mSample$tail_length,mSample$minu_nosite),
		"plus_wisite" = rep(mSample$tail_length,mSample$plus_wisite),
		"minu_wisite" = rep(mSample$tail_length,mSample$minu_wisite),
		"plus_nosite" = rep(mSample$tail_length,mSample$plus_nosite))


	pvals <- rbind(combn(4,2),apply(combn(4,2),2,function(x){
				ks.test(mDataSim[,x[1]],mDataSim[,x[2]], alternative = "greater")[[2]]
				}))
	print(pvals)
	
	print(round(apply(m[,-1],2,function(z){weighted.median(x = 0:250, w = z)})))
	print(args[3])
	print("fraction of mRNAs < 50 nt:")
	FracVec = apply(m[,-1],2,function(x){
		sum(x[8:21])/sum(x)
		})
	print(FracVec["plus_wisite"] - max(FracVec[names(FracVec) != "plus_wisite"]))
	m_melt<-melt(m,id="tail_length")
	
	cols <- c("black","grey","blue","red")
	m_melt$variable <- factor(m_melt$variable,levels = c("minu_nosite",
														"plus_nosite",
														"minu_wisite",
														"plus_wisite"))

	p1<-ggplot(m_melt,aes(x=tail_length,y=value,color=variable)) +
		geom_step(size = 0.5) +
		scale_x_continuous(name=NULL,expand=c(0,0), limits = c(0,50)) +
		scale_y_continuous(name=NULL,expand=c(0,0), limits = c(0,0.01)) +
		scale_colour_manual(values=cols) +
		theme_tim()


	
	ggsave(p1,file = paste0("figures/other/otherV9/HybridS1Histm_",args[3]),width=1.5,height=1.5,useDingbats=FALSE)}	
