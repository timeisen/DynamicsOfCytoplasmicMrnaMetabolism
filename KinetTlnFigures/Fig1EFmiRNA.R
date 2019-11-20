#Histograms of miRNA effects, miRNA fig 1EF
library(plyr)
library(reshape2)
library(tidyverse)
source("/lab/solexa_bartel/teisen/RNAseq/Scripts/general/ggplot_theme.R")
#Histograms of miRNAs.
##Uridylation

# args <- commandArgs(trailingOnly=TRUE)

#This code was modified on 2019 02 27 to create CDFs from direct ligation data. 
#It uses length matched no site data in the following folder: /lab/solexa_bartel/teisen/Tail-seq/PalSeqDirectLig/Analysis/miR_analysis/tpUTRlenMatching

for(i in 1:4){

	if(i == 1){
	PAL_minu_nosite = read_tsv("/lab/solexa_bartel/teisen/Tail-seq/PalSeqDirectLig/Analysis/miR_analysis_uridylation/miR1_uninduced_nosite_sitematchU.txt")
	PAL_plus_wisite = read_tsv("/lab/solexa_bartel/teisen/Tail-seq/PalSeqDirectLig/Analysis/miR_analysis_uridylation/miR1_induced_siteU.txt")
	PAL_minu_wisite = read_tsv("/lab/solexa_bartel/teisen/Tail-seq/PalSeqDirectLig/Analysis/miR_analysis_uridylation/miR1_uninduced_siteU.txt")
	PAL_plus_nosite = read_tsv("/lab/solexa_bartel/teisen/Tail-seq/PalSeqDirectLig/Analysis/miR_analysis_uridylation/miR1_induced_nosite_sitematchU.txt")
	}

	if(i == 2){
	PAL_minu_nosite = read_tsv("/lab/solexa_bartel/teisen/Tail-seq/PalSeqDirectLig/Analysis/miR_analysis_uridylation/miR1_uninduced_nosite_topsitematchU.txt")
	PAL_plus_wisite = read_tsv("/lab/solexa_bartel/teisen/Tail-seq/PalSeqDirectLig/Analysis/miR_analysis_uridylation/miR1_induced_topsiteU.txt")
	PAL_minu_wisite = read_tsv("/lab/solexa_bartel/teisen/Tail-seq/PalSeqDirectLig/Analysis/miR_analysis_uridylation/miR1_uninduced_topsiteU.txt")
	PAL_plus_nosite = read_tsv("/lab/solexa_bartel/teisen/Tail-seq/PalSeqDirectLig/Analysis/miR_analysis_uridylation/miR1_induced_nosite_topsitematchU.txt")
	}

	if(i == 3){
	PAL_minu_nosite = read_tsv("/lab/solexa_bartel/teisen/Tail-seq/PalSeqDirectLig/Analysis/miR_analysis_uridylation/miR155_uninduced_nosite_sitematchU.txt")
	PAL_plus_wisite = read_tsv("/lab/solexa_bartel/teisen/Tail-seq/PalSeqDirectLig/Analysis/miR_analysis_uridylation/miR155_induced_siteU.txt")
	PAL_minu_wisite = read_tsv("/lab/solexa_bartel/teisen/Tail-seq/PalSeqDirectLig/Analysis/miR_analysis_uridylation/miR155_uninduced_siteU.txt")
	PAL_plus_nosite = read_tsv("/lab/solexa_bartel/teisen/Tail-seq/PalSeqDirectLig/Analysis/miR_analysis_uridylation/miR155_induced_nosite_sitematchU.txt")
	}


	if(i == 4){
	PAL_minu_nosite = read_tsv("/lab/solexa_bartel/teisen/Tail-seq/PalSeqDirectLig/Analysis/miR_analysis_uridylation/miR155_uninduced_nosite_topsitematchU.txt")
	PAL_plus_wisite = read_tsv("/lab/solexa_bartel/teisen/Tail-seq/PalSeqDirectLig/Analysis/miR_analysis_uridylation/miR155_induced_topsiteU.txt")
	PAL_minu_wisite = read_tsv("/lab/solexa_bartel/teisen/Tail-seq/PalSeqDirectLig/Analysis/miR_analysis_uridylation/miR155_uninduced_topsiteU.txt")
	PAL_plus_nosite = read_tsv("/lab/solexa_bartel/teisen/Tail-seq/PalSeqDirectLig/Analysis/miR_analysis_uridylation/miR155_induced_nosite_topsitematchU.txt")
	}
	


	PAL_minu_nosite <- mutate(PAL_minu_nosite, minu_nosite = rowSums(PAL_minu_nosite[,2:11])/rowSums(PAL_minu_nosite[,-1])) %>%
		select(tail_length,minu_nosite) %>% 
		arrange(tail_length)
	
	PAL_plus_wisite <- mutate(PAL_plus_wisite, plus_wisite = rowSums(PAL_plus_wisite[,2:11])/rowSums(PAL_plus_wisite[,-1])) %>%
		select(tail_length,plus_wisite) %>% 
		arrange(tail_length)
		
	PAL_minu_wisite <- mutate(PAL_minu_wisite, minu_wisite = rowSums(PAL_minu_wisite[,2:11])/rowSums(PAL_minu_wisite[,-1])) %>%
		select(tail_length,minu_wisite) %>% 
		arrange(tail_length)
		
	PAL_plus_nosite <- mutate(PAL_plus_nosite, plus_nosite = rowSums(PAL_plus_nosite[,2:11])/rowSums(PAL_plus_nosite[,-1])) %>%
		select(tail_length,plus_nosite) %>% 
		arrange(tail_length)
	
	#Normalization of the rows that have values above the median value. 
	# PAL_minu_nosite_median<-median(rowSums(PAL_minu_nosite[-1]))
	# PAL_plus_wisite_median<-median(rowSums(PAL_plus_wisite[-1]))
	# PAL_minu_wisite_median<-median(rowSums(PAL_minu_wisite[-1]))
	# PAL_plus_nosite_median<-median(rowSums(PAL_plus_nosite[-1]))
	
	
	# PAL_minu_nosite_high<-which(rowSums(PAL_minu_nosite[-1])>PAL_minu_nosite_median)
	# PAL_plus_wisite_high<-which(rowSums(PAL_plus_wisite[-1])>PAL_plus_wisite_median)
	# PAL_minu_wisite_high<-which(rowSums(PAL_minu_wisite[-1])>PAL_minu_wisite_median)
	# PAL_plus_nosite_high<-which(rowSums(PAL_plus_nosite[-1])>PAL_plus_nosite_median)
	
	# PAL_minu_nosite[PAL_minu_nosite_high,-1] = PAL_minu_nosite[PAL_minu_nosite_high,-1]/rowSums(PAL_minu_nosite[PAL_minu_nosite_high,-1])*	PAL_minu_nosite_median
	# PAL_plus_wisite[PAL_plus_wisite_high,-1] = PAL_plus_wisite[PAL_plus_wisite_high,-1]/rowSums(PAL_plus_wisite[PAL_plus_wisite_high,-1])*	PAL_plus_wisite_median
	# PAL_minu_wisite[PAL_minu_wisite_high,-1] = PAL_minu_wisite[PAL_minu_wisite_high,-1]/rowSums(PAL_minu_wisite[PAL_minu_wisite_high,-1])*	PAL_minu_wisite_median
	# PAL_plus_nosite[PAL_plus_nosite_high,-1] = PAL_plus_nosite[PAL_plus_nosite_high,-1]/rowSums(PAL_plus_nosite[PAL_plus_nosite_high,-1])*	PAL_plus_nosite_median
	
	m <- reduce(list(PAL_minu_nosite, 
		PAL_plus_wisite, 
		PAL_minu_wisite, 
		PAL_plus_nosite), inner_join, by = "tail_length")
	print(m)
	m_melt<-melt(m,id="tail_length")
	
	cols <- c("black","grey","blue","red")
	m_melt$variable <- factor(m_melt$variable,levels = c("minu_nosite",
														"plus_nosite",
														"minu_wisite",
														"plus_wisite"))

	p1<-ggplot(m_melt,aes(x=tail_length,y=value,color=variable)) +
		geom_step(size = 0.5) +
		scale_x_continuous(name=NULL,expand=c(0,0), limits = c(0,50)) +
		scale_y_continuous(name=NULL,expand=c(0,0), limits = c(0,1)) +
		scale_colour_manual(values=cols) +
		theme_tim()


	
	ggsave(p1,file = paste0("figures/other/otherV9/HybridS1Histm_",i,".pdf"),width=1.5,height=1.5,useDingbats=FALSE)}	
