library(plyr)
library(reshape2)
library(tidyverse)
#library(data.table)

hl<-2.249019

pl<-NULL
pul<-NULL

# pul<-1/((1-exp(c(2/3,1,2,4,8)*log(2)/-2.24))/(1-exp(1*log(2)/-2.24)))*.1

pul<-rep(0.003,5)

pl[1]=1
pl[2]=1
pl[3]=1
pl[4]=1
pl[5]=1

#masses for normalization
masses = c(
	500,
	500,
	250,
	200,
	100,
	50)


#Add SS
stds<-read.table("/lab/solexa_bartel/teisen/RNAseq/Annotation_files/TAIL_seq/list_of_standards_SWE_notation.txt")
stds$V1<-as.character(stds$V1)
quant_stds<-stds$V1[1:11]

det_SS_norm<-function(tp8hr,tpSS_,tp8hr_norm_counts,tpSS__norm_counts){
	hl<-read.table("processed_files/halflifeComparisons/miR-155_minu_polyA_halflives_logspace_global_offset_noSS.txt",head = TRUE)
	hl<-hl[,c(1,2,5)]
	colnames(hl)<-c("accession","tx_rate","halflife")
	short_hl<-as.character(hl<-hl[which(hl$halflife < 40/60 & hl$halflife > 20/60),]$accession)
	tp8hr$rs<-rowSums(tp8hr[,-1]/sum(tp8hr[,-1])/tp8hr_norm_counts)
	tpSS_$rs<-rowSums(tpSS_[,-1]/sum(tpSS_[,-1])/tpSS__norm_counts)
	tp8hr<-tp8hr[,c(1,253)]
	tpSS_<-tpSS_[,c(1,253)]
	
	m<-merge(tp8hr,tpSS_,by="accession")
	colnames(m)[c(2,3)]<-c("tp8hr","tpSS_")
	m<-m[which(m$accession %in% short_hl),]
	return(median(m[,2])/median(m[,3]))}

ShortTailNorm <- function(rowX){ #For figuring out which genes have no reads in their first 8 nt to prevent NAs. 
	if (rowX["fracShort"] != 0) {rowX[1:8] = rowX[253:260]/sum(rowX[253:260])*rowX["fracShort"]} #one shifted left because no accession code
	else {rowX[1:8] = 0}
	return(rowX)
}

for(i in 1:10){
	print("Read files in ...")
	tp40m<-read.table(paste0("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/single_tag_files/bootstrap/40min_minus_tag_reformat_no_cutoff_miR-155_BOOTSTRAP_",i,".txt"),sep="\t")
	tp1hr<-read.table(paste0("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/single_tag_files/bootstrap/1hr_minus_tag_reformat_no_cutoff_miR-155_BOOTSTRAP_",i,".txt"),sep="\t")
	tp2hr<-read.table(paste0("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/single_tag_files/bootstrap/2hr_minus_tag_reformat_no_cutoff_miR-155_BOOTSTRAP_",i,".txt"),sep="\t")
	tp4hr<-read.table(paste0("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/single_tag_files/bootstrap/4hr_minus_tag_reformat_no_cutoff_miR-155_BOOTSTRAP_",i,".txt"),sep="\t")
	tp8hr<-read.table(paste0("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/single_tag_files/bootstrap/8hr_minus_tag_reformat_no_cutoff_miR-155_BOOTSTRAP_",i,".txt"),sep="\t")
	tpSS_<-read.table(paste0("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/single_tag_files/bootstrap/SS_minus_tag_reformat_no_cutoff_miR-155_BOOTSTRAP_",i,".txt"),sep="\t")
	tpSWE<-read.table(paste0("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/single_tag_files/bootstrap/tail_lengths_annot_end_only_CLEANED20190731_reformat_BOOTSTRAP_",i,".txt"),sep="\t")

	# tpSWE<-read.table(paste0("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/single_tag_files/bootstrap/tail_lengths_annot_end_only_CLEANED20190619_reformat_BOOTSTRAP_",i,".txt"),sep="\t")

	print(nrow(tp40m[which(rowSums(tp40m[,-1])>=50),]))
	print(nrow(tp1hr[which(rowSums(tp1hr[,-1])>=50),]))
	print(nrow(tp2hr[which(rowSums(tp2hr[,-1])>=50),]))
	print(nrow(tp4hr[which(rowSums(tp4hr[,-1])>=50),]))
	print(nrow(tp8hr[which(rowSums(tp8hr[,-1])>=50),]))
	print(nrow(tpSS_[which(rowSums(tpSS_[,-1])>=50),]))
	print(nrow(tpSWE[which(rowSums(tpSWE[,-1])>=50),]))

	tpSWE <- tpSWE[which(rowSums(tpSWE[,-1])>=50),]

	cutoff_tp40m <- tp40m[which(rowSums(tp40m[,-1])==50),1]
	cutoff_tp1hr <- tp1hr[which(rowSums(tp1hr[,-1])==50),1] #sometimes a gene is not in the steady state dataset by chance
	cutoff_tp2hr <- tp2hr[which(rowSums(tp2hr[,-1])==50),1]
	cutoff_tp4hr <- tp4hr[which(rowSums(tp4hr[,-1])==50),1]
	cutoff_tp8hr <- tp8hr[which(rowSums(tp8hr[,-1])==50),1]
	cutoff_tpSS_ <- tpSS_[which(rowSums(tpSS_[,-1])==50),1]

	print(length(cutoff_tp40m))
	print(length(cutoff_tp1hr))
	print(length(cutoff_tp2hr))
	print(length(cutoff_tp4hr))
	print(length(cutoff_tp8hr))
	print(length(cutoff_tpSS_))

	colnames(tp40m)[1] <- "accession"
	colnames(tp1hr)[1] <- "accession"
	colnames(tp2hr)[1] <- "accession"
	colnames(tp4hr)[1] <- "accession"
	colnames(tp8hr)[1] <- "accession"
	colnames(tpSS_)[1] <- "accession"

	tp40m$accession<-as.character(tp40m$accession)
	tp1hr$accession<-as.character(tp1hr$accession)
	tp2hr$accession<-as.character(tp2hr$accession)
	tp4hr$accession<-as.character(tp4hr$accession)
	tp8hr$accession<-as.character(tp8hr$accession)
	tpSS_$accession<-as.character(tpSS_$accession)



	#any gene that has fewer than 50 tags in the ss or tp datasets is removed later
	# tpSS__genes_below_cutoff<-tpSS_[which(rowSums(tpSS_[-1])<50),1]
	# tp40m_genes_below_cutoff<-c(tp40m[which(rowSums(tp40m[-1])<50),1],tpSS__genes_below_cutoff)
	# tp1hr_genes_below_cutoff<-c(tp1hr[which(rowSums(tp1hr[-1])<50),1],tpSS__genes_below_cutoff)
	# tp2hr_genes_below_cutoff<-c(tp2hr[which(rowSums(tp2hr[-1])<50),1],tpSS__genes_below_cutoff)
	# tp4hr_genes_below_cutoff<-c(tp4hr[which(rowSums(tp4hr[-1])<50),1],tpSS__genes_below_cutoff)
	# tp8hr_genes_below_cutoff<-c(tp8hr[which(rowSums(tp8hr[-1])<50),1],tpSS__genes_below_cutoff)

	tp40m_norm_counts<-sum(rowSums(tp40m[which(tp40m$accession %in% quant_stds),-1]))/sum(rowSums(tp40m[which(!tp40m$accession %in% stds$V1),-1]))*masses[1]
	tp1hr_norm_counts<-sum(rowSums(tp1hr[which(tp1hr$accession %in% quant_stds),-1]))/sum(rowSums(tp1hr[which(!tp1hr$accession %in% stds$V1),-1]))*masses[2]
	tp2hr_norm_counts<-sum(rowSums(tp2hr[which(tp2hr$accession %in% quant_stds),-1]))/sum(rowSums(tp2hr[which(!tp2hr$accession %in% stds$V1),-1]))*masses[3]
	tp4hr_norm_counts<-sum(rowSums(tp4hr[which(tp4hr$accession %in% quant_stds),-1]))/sum(rowSums(tp4hr[which(!tp4hr$accession %in% stds$V1),-1]))*masses[4]
	tp8hr_norm_counts<-sum(rowSums(tp8hr[which(tp8hr$accession %in% quant_stds),-1]))/sum(rowSums(tp8hr[which(!tp8hr$accession %in% stds$V1),-1]))*masses[5]
	tpSS__norm_counts<-sum(rowSums(tpSS_[which(tpSS_$accession %in% quant_stds),-1]))/sum(rowSums(tpSS_[which(!tpSS_$accession %in% stds$V1),-1]))*masses[6]
	print(tpSS__norm_counts)
	tpSS__norm_counts<-tpSS__norm_counts/det_SS_norm(tp8hr,tpSS_,tp8hr_norm_counts,tpSS__norm_counts)
	print(tpSS__norm_counts)

	#remove stds here?
	tp40m <- tp40m[which(!tp40m$accession %in% stds$V1),]
	tp1hr <- tp1hr[which(!tp1hr$accession %in% stds$V1),]
	tp2hr <- tp2hr[which(!tp2hr$accession %in% stds$V1),]
	tp4hr <- tp4hr[which(!tp4hr$accession %in% stds$V1),]
	tp8hr <- tp8hr[which(!tp8hr$accession %in% stds$V1),]
	tpSS_ <- tpSS_[which(!tpSS_$accession %in% stds$V1),]

	# tpSS_<-tpSS_[which(rowSums(tpSS_[,-1])>=50),]
	tp40m_merge<-merge(tp40m,tpSS_,by="accession")
	tp1hr_merge<-merge(tp1hr,tpSS_,by="accession")
	tp2hr_merge<-merge(tp2hr,tpSS_,by="accession")
	tp4hr_merge<-merge(tp4hr,tpSS_,by="accession")
	tp8hr_merge<-merge(tp8hr,tpSS_,by="accession")

	print("Subtraction ...")
	tp40m_merge[,2:252] <- tp40m_merge[,2:252]/sum(tp40m_merge[,2:252])/tp40m_norm_counts
	tp1hr_merge[,2:252] <- tp1hr_merge[,2:252]/sum(tp1hr_merge[,2:252])/tp1hr_norm_counts
	tp2hr_merge[,2:252] <- tp2hr_merge[,2:252]/sum(tp2hr_merge[,2:252])/tp2hr_norm_counts
	tp4hr_merge[,2:252] <- tp4hr_merge[,2:252]/sum(tp4hr_merge[,2:252])/tp4hr_norm_counts
	tp8hr_merge[,2:252] <- tp8hr_merge[,2:252]/sum(tp8hr_merge[,2:252])/tp8hr_norm_counts

	tp40m_merge[,253:503] <- tp40m_merge[,253:503]/sum(tp40m_merge[,253:503])/tpSS__norm_counts
	tp1hr_merge[,253:503] <- tp1hr_merge[,253:503]/sum(tp1hr_merge[,253:503])/tpSS__norm_counts
	tp2hr_merge[,253:503] <- tp2hr_merge[,253:503]/sum(tp2hr_merge[,253:503])/tpSS__norm_counts
	tp4hr_merge[,253:503] <- tp4hr_merge[,253:503]/sum(tp4hr_merge[,253:503])/tpSS__norm_counts
	tp8hr_merge[,253:503] <- tp8hr_merge[,253:503]/sum(tp8hr_merge[,253:503])/tpSS__norm_counts

	tp40m_bs <- (tp40m_merge[,253:503]*pul[1]-tp40m_merge[,2:252])/(pul[1]-pl[1])
	tp1hr_bs <- (tp1hr_merge[,253:503]*pul[2]-tp1hr_merge[,2:252])/(pul[2]-pl[2])
	tp2hr_bs <- (tp2hr_merge[,253:503]*pul[3]-tp2hr_merge[,2:252])/(pul[3]-pl[3])
	tp4hr_bs <- (tp4hr_merge[,253:503]*pul[4]-tp4hr_merge[,2:252])/(pul[4]-pl[4])
	tp8hr_bs <- (tp8hr_merge[,253:503]*pul[5]-tp8hr_merge[,2:252])/(pul[5]-pl[5])
	tpSS__bs <- tpSS_[-1]/sum(tpSS_[-1])/tpSS__norm_counts

	print("Applying cutoffs ...")

	tp40m_bs$accession <- tp40m_merge$accession
	tp1hr_bs$accession <- tp1hr_merge$accession
	tp2hr_bs$accession <- tp2hr_merge$accession
	tp4hr_bs$accession <- tp4hr_merge$accession
	tp8hr_bs$accession <- tp8hr_merge$accession
	tpSS__bs$accession <- tpSS_$accession

	print("Number of genes total:")
	print(nrow(tp40m_bs))
	print(nrow(tp1hr_bs))
	print(nrow(tp2hr_bs))
	print(nrow(tp4hr_bs))
	print(nrow(tp8hr_bs))
	print(nrow(tpSS__bs))

	cutoff_tp40m <- mean(rowSums(tp40m_bs[which(tp40m_bs$accession %in% cutoff_tp40m),-252]))
	cutoff_tp1hr <- mean(rowSums(tp1hr_bs[which(tp1hr_bs$accession %in% cutoff_tp1hr),-252]))
	cutoff_tp2hr <- mean(rowSums(tp2hr_bs[which(tp2hr_bs$accession %in% cutoff_tp2hr),-252]))
	cutoff_tp4hr <- mean(rowSums(tp4hr_bs[which(tp4hr_bs$accession %in% cutoff_tp4hr),-252]))
	cutoff_tp8hr <- mean(rowSums(tp8hr_bs[which(tp8hr_bs$accession %in% cutoff_tp8hr),-252]))
	cutoff_tpSS_ <- mean(rowSums(tpSS__bs[which(tpSS__bs$accession %in% cutoff_tpSS_),-252]))

	tp40m_bs<-tp40m_bs[which(rowSums(tp40m_bs[,-252])>=cutoff_tp40m),]
	tp1hr_bs<-tp1hr_bs[which(rowSums(tp1hr_bs[,-252])>=cutoff_tp1hr),]
	tp2hr_bs<-tp2hr_bs[which(rowSums(tp2hr_bs[,-252])>=cutoff_tp2hr),]
	tp4hr_bs<-tp4hr_bs[which(rowSums(tp4hr_bs[,-252])>=cutoff_tp4hr),]
	tp8hr_bs<-tp8hr_bs[which(rowSums(tp8hr_bs[,-252])>=cutoff_tp8hr),]
	tpSS__bs<-tpSS__bs[which(rowSums(tpSS__bs[,-252])>=cutoff_tpSS_),]

	print("Number of genes above cutoffs:")
	print(nrow(tp40m_bs))
	print(nrow(tp1hr_bs))
	print(nrow(tp2hr_bs))
	print(nrow(tp4hr_bs))
	print(nrow(tp8hr_bs))
	print(nrow(tpSS__bs))

	##########################################################################################
	#Lines added on 2019 04 11 to rescale the data based on the RNAseq libraries. 
	
	###Note that it doesn't matter what dataset is used here as long as the time points (not the steady state) 
	#are the based on the RNAseq scaling. The steady state will be rescaled later.
	expr_data <- read.table("processed_files/rnaseqDatasets/miR-155_minu_polyA_30minHLnorm.txt",head = TRUE,sep = "\t")
	expr_norm <- colSums(expr_data[,-1])
	###
	
	
	#This is complicated and needs to be cleaned up. 
	
	data_total_sum <- sum(
		tp40m_bs[,-252],
		tp1hr_bs[,-252],
		tp2hr_bs[,-252],
		tp4hr_bs[,-252],
		tp8hr_bs[,-252],
		tpSS__bs[,-252]
		)
	
	tp40m_bs[,-252] <- tp40m_bs[,-252]/sum(tp40m_bs[,-252]) * expr_norm[1] 
	tp1hr_bs[,-252] <- tp1hr_bs[,-252]/sum(tp1hr_bs[,-252]) * expr_norm[2] 
	tp2hr_bs[,-252] <- tp2hr_bs[,-252]/sum(tp2hr_bs[,-252]) * expr_norm[3] 
	tp4hr_bs[,-252] <- tp4hr_bs[,-252]/sum(tp4hr_bs[,-252]) * expr_norm[4] 
	tp8hr_bs[,-252] <- tp8hr_bs[,-252]/sum(tp8hr_bs[,-252]) * expr_norm[5] 
	tpSS__bs[,-252] <- tpSS__bs[,-252]/sum(tpSS__bs[,-252]) * expr_norm[6] 
	
	data_total_sum_prenorm <- sum(
		tp40m_bs[,-252],
		tp1hr_bs[,-252],
		tp2hr_bs[,-252],
		tp4hr_bs[,-252],
		tp8hr_bs[,-252],
		tpSS__bs[,-252]
		)
	
	tp40m_bs[,-252] <- tp40m_bs[,-252] / data_total_sum_prenorm * data_total_sum
	tp1hr_bs[,-252] <- tp1hr_bs[,-252] / data_total_sum_prenorm * data_total_sum
	tp2hr_bs[,-252] <- tp2hr_bs[,-252] / data_total_sum_prenorm * data_total_sum
	tp4hr_bs[,-252] <- tp4hr_bs[,-252] / data_total_sum_prenorm * data_total_sum
	tp8hr_bs[,-252] <- tp8hr_bs[,-252] / data_total_sum_prenorm * data_total_sum
	tpSS__bs[,-252] <- tpSS__bs[,-252] / data_total_sum_prenorm * data_total_sum
	
	##########################################################################################

	# tp40m_bs<-tp40m_bs[which(!tp40m_bs$accession %in% stds$V1),]
	# tp1hr_bs<-tp1hr_bs[which(!tp1hr_bs$accession %in% stds$V1),]
	# tp2hr_bs<-tp2hr_bs[which(!tp2hr_bs$accession %in% stds$V1),]
	# tp4hr_bs<-tp4hr_bs[which(!tp4hr_bs$accession %in% stds$V1),]
	# tp8hr_bs<-tp8hr_bs[which(!tp8hr_bs$accession %in% stds$V1),]
	# tpSS__bs<-tpSS__bs[which(!tpSS__bs$accession %in% stds$V1),]
	#MODIFY THE STEADY STATE DATA TO INCORPORATE TAIL-SEQ

	#PREVIOUS CODE
	# tpSWE$rs = rowSums(tpSWE[,-1])
	# mSWEaccession = as.character(tpSWE[,1])
	# tpSWE = tpSWE[,2:9]/tpSWE$rs
	# tpSWE$accession = mSWEaccession
	# tpSS__bs$accession = as.character(tpSS__bs$accession)
	# mPALSWE <- merge(tpSS__bs,tpSWE,by="accession")
	# # break #these lines need to be fixed
	# mPALSWE[,2:9] = rowSums(mPALSWE[,2:252])*mPALSWE[,253:260]
	# mPALSWE[,253:260] <- NULL
	# tpSS__bs <- mPALSWE[,c(2:252,1)]

	DIR <- tpSWE
 	SPL <- tpSS__bs
 	colnames(SPL) <- c(paste0("SPL_",0:250),"accession")
 	colnames(DIR) <- c("accession",paste0("DIR_",0:250))
 	SPL <- SPL[,c(252,1:251)] #make accession the first column

	DIR$fracShort   = rowSums(DIR[,2:9])/rowSums(DIR[,2:252])
	DIR$fracLong    = rowSums(DIR[,10:252])/rowSums(DIR[,2:252])
	SPL$expr        = rowSums(SPL[,-1])
	SPLDIR          = inner_join(SPL,DIR,by="accession")
	SPLDIR[,-1]     = t(apply(SPLDIR[,-1],1,ShortTailNorm)) #Short tailed norm, revised. 
	# SPLDIR[,2:9]    = SPLDIR[,254:261]/rowSums(SPLDIR[,254:261])*SPLDIR$fracShort #short tail normalization
	SPLDIR[,10:252] = SPLDIR[,10:252]/rowSums(SPLDIR[,10:252])*SPLDIR$fracLong #long tail normalization
	SPLDIR[,2:252]  = SPLDIR[,2:252]/rowSums(SPLDIR[,2:252])*SPLDIR$expr #expression normalization
	
	tpSS__bs <- SPLDIR[,c(2:252,1)] #switch the accession code back to the last column
 
	print("Writing files ...")
	write.table(tp40m_bs,file=paste0("single_tag_files/background_subtracted_single_tag_files/bootstrap/tp40m_bs_minus_miR-155_v7_BOOTSTRAP_HYBRID20190731_",i,".txt"),quote=FALSE,row.names=FALSE,sep="\t",col.names=FALSE)
	write.table(tp1hr_bs,file=paste0("single_tag_files/background_subtracted_single_tag_files/bootstrap/tp1hr_bs_minus_miR-155_v7_BOOTSTRAP_HYBRID20190731_",i,".txt"),quote=FALSE,row.names=FALSE,sep="\t",col.names=FALSE)
	write.table(tp2hr_bs,file=paste0("single_tag_files/background_subtracted_single_tag_files/bootstrap/tp2hr_bs_minus_miR-155_v7_BOOTSTRAP_HYBRID20190731_",i,".txt"),quote=FALSE,row.names=FALSE,sep="\t",col.names=FALSE)
	write.table(tp4hr_bs,file=paste0("single_tag_files/background_subtracted_single_tag_files/bootstrap/tp4hr_bs_minus_miR-155_v7_BOOTSTRAP_HYBRID20190731_",i,".txt"),quote=FALSE,row.names=FALSE,sep="\t",col.names=FALSE)
	write.table(tp8hr_bs,file=paste0("single_tag_files/background_subtracted_single_tag_files/bootstrap/tp8hr_bs_minus_miR-155_v7_BOOTSTRAP_HYBRID20190731_",i,".txt"),quote=FALSE,row.names=FALSE,sep="\t",col.names=FALSE)
	write.table(tpSS__bs,file=paste0("single_tag_files/background_subtracted_single_tag_files/bootstrap/tpSS_minus_norm_miR-155_v7_BOOTSTRAP_HYBRID20190731_",i,".txt"),quote=FALSE,row.names=FALSE,sep="\t",col.names=FALSE)

	#_____________________________
	#_____________________________
	print("Writing halflife fitting files ...")
	tp40m_bs_rs <- data.frame(tp40m_bs[,252],rowSums(tp40m_bs[,-c(252)]))
	tp1hr_bs_rs <- data.frame(tp1hr_bs[,252],rowSums(tp1hr_bs[,-c(252)]))
	tp2hr_bs_rs <- data.frame(tp2hr_bs[,252],rowSums(tp2hr_bs[,-c(252)]))
	tp4hr_bs_rs <- data.frame(tp4hr_bs[,252],rowSums(tp4hr_bs[,-c(252)]))
	tp8hr_bs_rs <- data.frame(tp8hr_bs[,252],rowSums(tp8hr_bs[,-c(252)]))
	tpSS__bs_rs <- data.frame(tpSS__bs[,252],rowSums(tpSS__bs[,-c(252)]))
	
	colnames(tp40m_bs_rs) <- c("accession","x40m")
	colnames(tp1hr_bs_rs) <- c("accession","x1hr")
	colnames(tp2hr_bs_rs) <- c("accession","x2hr")
	colnames(tp4hr_bs_rs) <- c("accession","x4hr")
	colnames(tp8hr_bs_rs) <- c("accession","x8hr")
	colnames(tpSS__bs_rs) <- c("accession","xSS_")
	
	mRS<-join_all(list(
		tp40m_bs_rs,
		tp1hr_bs_rs,
		tp2hr_bs_rs,
		tp4hr_bs_rs,
		tp8hr_bs_rs,
		tpSS__bs_rs),by="accession",type="full")
	
	mRS <- mRS[complete.cases(mRS),]
	mRS <- mRS[which(!mRS$x40m < 0),] #no negative values
	
	print(nrow(mRS))
	write.table(mRS,file = paste0("processed_files/palSeqExprNorm/bootstrap/miR-155_minu_PAL_expr_norm_RNAseq_scaled_v7_20190731_bootstrap_",i,".txt"),quote = FALSE, sep = "\t",row.names = FALSE)
	
	#_____________________________
	#_____________________________

	print("Fraction of negative values:")
	print(length(which(tp40m_bs[,-252]<0))/nrow(tp40m_bs[,-252])/ncol(tp40m_bs[,-252]))
	print(length(which(tp1hr_bs[,-252]<0))/nrow(tp1hr_bs[,-252])/ncol(tp1hr_bs[,-252]))
	print(length(which(tp2hr_bs[,-252]<0))/nrow(tp2hr_bs[,-252])/ncol(tp2hr_bs[,-252]))
	print(length(which(tp4hr_bs[,-252]<0))/nrow(tp4hr_bs[,-252])/ncol(tp4hr_bs[,-252]))
	print(length(which(tp8hr_bs[,-252]<0))/nrow(tp8hr_bs[,-252])/ncol(tp8hr_bs[,-252]))
	print(length(which(tpSS__bs[,-252]<0))/nrow(tpSS__bs[,-252])/ncol(tpSS__bs[,-252]))


	print("Scaled abundances:")
	print(sum(tp40m_bs[,-252]))
	print(sum(tp1hr_bs[,-252]))
	print(sum(tp2hr_bs[,-252]))
	print(sum(tp4hr_bs[,-252]))
	print(sum(tp8hr_bs[,-252]))
	print(sum(tpSS__bs[,-252]))


	v<-c(
	sum(tp40m_bs[,-252]),
	sum(tp1hr_bs[,-252]),
	sum(tp2hr_bs[,-252]),
	sum(tp4hr_bs[,-252]),
	sum(tp8hr_bs[,-252]),
	sum(tpSS__bs[,-252])
	)

	probs = c(0.01,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.99)

	tp40m_m<-data.frame(apply(tp40m_bs[,-252],1,function(x){weighted.mean(0:250,w=x)}))
	tp1hr_m<-data.frame(apply(tp1hr_bs[,-252],1,function(x){weighted.mean(0:250,w=x)}))
	tp2hr_m<-data.frame(apply(tp2hr_bs[,-252],1,function(x){weighted.mean(0:250,w=x)}))
	tp4hr_m<-data.frame(apply(tp4hr_bs[,-252],1,function(x){weighted.mean(0:250,w=x)}))
	tp8hr_m<-data.frame(apply(tp8hr_bs[,-252],1,function(x){weighted.mean(0:250,w=x)}))
	tpSS__m<-data.frame(apply(tpSS__bs[,-252],1,function(x){weighted.mean(0:250,w=x)}))


	tp40m_m$accession = tp40m_bs[,252]
	tp1hr_m$accession = tp1hr_bs[,252]
	tp2hr_m$accession = tp2hr_bs[,252]
	tp4hr_m$accession = tp4hr_bs[,252]
	tp8hr_m$accession = tp8hr_bs[,252]
	tpSS__m$accession = tpSS__bs[,252]



	tp40m_m<-tp40m_m[,c(2,1)]
	tp1hr_m<-tp1hr_m[,c(2,1)]
	tp2hr_m<-tp2hr_m[,c(2,1)]
	tp4hr_m<-tp4hr_m[,c(2,1)]
	tp8hr_m<-tp8hr_m[,c(2,1)]
	tpSS__m<-tpSS__m[,c(2,1)]


	m<-join_all(list(
	tp40m_m,
	tp1hr_m,
	tp2hr_m,
	tp4hr_m,
	tp8hr_m,
	tpSS__m),by="accession",type="full")

	m<-m[complete.cases(m),]
	colnames(m)<-c("accession","x40m","x1hr","x2hr","x4hr","x8hr","xtpSS_")
	print("Means:")
	print(apply(m[,-1],2,function(x){quantile(x,probs=probs)}))
	# print("Writing mean files...")
	# write.table(m,"/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/processed_files/miR-155_minus_sample_mean_tails_50tags_background_subtracted_v5.txt",quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")

	expr<-read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/model_input_files/miR-155_minus_sample_mean_tails_50tags_with_expr.txt",head=FALSE)
	expr<-expr[,c(1,8:13)]

	colnames(expr)[1]<-"accession"

	m<-merge(m,expr)
	print("Correlations:")
	print(cor(m[-1],method="spearman"))

	# write.table(m,"/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/model_input_files/miR-155_minus_sample_mean_tails_50tags_with_expr_background_subtracted_v5.txt",quote=FALSE,row.names=FALSE,col.names=FALSE,sep="\t")
}