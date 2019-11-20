library(plyr)
library(ggplot2)
library(data.table)
library(tidyverse)

for(i in 1:10){
	mTemp <- read_tsv(paste0("rate_constant_measurements/bootstrap/miR-155_minus_miR-155_minus_samples_UNLINKV3_V81_H5_Run3_BOOTSTRAP_",i,".txt"),col_names = FALSE)
	colnames(mTemp)[c(1,6)] = c("accession","residual")

	bootstrap_params <- read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/rate_constant_measurements/global_model_fits/bootstrapHYBRID/miR-155_minus_rand_g_rates_plogis_unlinked_HYBRID20190801_compiled_summarize.txt",sep="\t",head=TRUE)

	pLogEl = plogis(230, loc = bootstrap_params[i,]$loc,scale = bootstrap_params[i,]$sc)


	m155 <- data.frame(mTemp[order(mTemp$residual),] %>% distinct(accession,.keep_all = TRUE)) #keep the lowest residual.
	colnames(m155)<-c("accession","st_","a_","k_","b_","r_")

	#b
	m155[which(m155$b_*pLogEl < 0.003),]$b_ = 0.003/pLogEl
	m155[which(m155$b_*pLogEl > 3),]$b_ = 3/pLogEl

	#a
	m155[which(m155$a_ < 1E-8),]$a_ = 1E-8
	m155[which(m155$a_ > 1E-5),]$a_ = 1E-5

	#k
	m155[which(m155$k_ < 0.1),]$k_ = 0.1
	m155[which(m155$k_ > 100),]$k_ = 100

	write_tsv(m155, path = paste0("rate_constant_measurements/bootstrap/miR-155_minus_miR-155_minus_samples_UNLINKV3_V81_H5_Run3_BOOTSTRAP_",i,"_TRUNC.txt"))

}
