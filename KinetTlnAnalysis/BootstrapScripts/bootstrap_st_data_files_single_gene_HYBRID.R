#This script generates the single tails for model inputs
library(plyr)

for(i in 1:10){
	#keep steady state in there?
	tp40m_ <- read.table(paste0("single_tag_files/background_subtracted_single_tag_files/bootstrap/tp40m_bs_minus_miR-155_v7_BOOTSTRAP_HYBRID20190731_",i,".txt"))
	tp1hr_ <- read.table(paste0("single_tag_files/background_subtracted_single_tag_files/bootstrap/tp1hr_bs_minus_miR-155_v7_BOOTSTRAP_HYBRID20190731_",i,".txt"))
	tp2hr_ <- read.table(paste0("single_tag_files/background_subtracted_single_tag_files/bootstrap/tp2hr_bs_minus_miR-155_v7_BOOTSTRAP_HYBRID20190731_",i,".txt"))
	tp4hr_ <- read.table(paste0("single_tag_files/background_subtracted_single_tag_files/bootstrap/tp4hr_bs_minus_miR-155_v7_BOOTSTRAP_HYBRID20190731_",i,".txt"))
	tp8hr_ <- read.table(paste0("single_tag_files/background_subtracted_single_tag_files/bootstrap/tp8hr_bs_minus_miR-155_v7_BOOTSTRAP_HYBRID20190731_",i,".txt"))
	SS____ <- read.table(paste0("single_tag_files/background_subtracted_single_tag_files/bootstrap/tpSS_minus_norm_miR-155_v7_BOOTSTRAP_HYBRID20190731_RescaleHybrid_",i,".txt"))
	# SS____[,253] <- SS____[,1]
	# SS____[,1] <- NULL
	
	#removes the 250 and the first 19 nt. 
	tp40m_<-tp40m_[,-c(1:20,251)]
	tp1hr_<-tp1hr_[,-c(1:20,251)]
	tp2hr_<-tp2hr_[,-c(1:20,251)]
	tp4hr_<-tp4hr_[,-c(1:20,251)]
	tp8hr_<-tp8hr_[,-c(1:20,251)]
	SS____<-SS____[,-c(251)]
	
	colnames(tp40m_)<-c(paste0("tp40m_",20:249),"accession")
	colnames(tp1hr_)<-c(paste0("tp1hr_",20:249),"accession")
	colnames(tp2hr_)<-c(paste0("tp2hr_",20:249),"accession")
	colnames(tp4hr_)<-c(paste0("tp4hr_",20:249),"accession")
	colnames(tp8hr_)<-c(paste0("tp8hr_",20:249),"accession")
	colnames(SS____)<-c(paste0("SS____",0:249),"accession")
	
	#scaling the datasets should already be scaled
	# tp40m_[,-251]<-tp40m_[,-251]/1
	# tp1hr_[,-251]<-tp1hr_[,-251]/1.74449207
	# tp2hr_[,-251]<-tp2hr_[,-251]/0.193279838
	# tp4hr_[,-251]<-tp4hr_[,-251]/0.12829349
	# tp8hr_[,-251]<-tp8hr_[,-251]/0.079662677
	# SS____[,-251]<-SS____[,-251]/0.048435329
	
	tp40m_ <- tp40m_[,231:1] #Modified to change it to 20. 
	tp1hr_ <- tp1hr_[,231:1] #Modified to change it to 20. 
	tp2hr_ <- tp2hr_[,231:1] #Modified to change it to 20. 
	tp4hr_ <- tp4hr_[,231:1] #Modified to change it to 20. 
	tp8hr_ <- tp8hr_[,231:1] #Modified to change it to 20. 
	SS____ <- SS____[,251:1]
	
	m<-join_all(list(
		tp40m_,
		tp1hr_,
		tp2hr_,
		tp4hr_,
		tp8hr_,
		SS____),by="accession",type="full")
	
	m<-m[complete.cases(m),]
	print(nrow(m))

	write.table(m,paste0("model_input_files/bootstrap/miR-155_minus_sample_mean_tails_50tags_with_expr_background_subtracted_v7_st_BOOTSTRAP_HYBRID20190801_Rescale_",i,".txt"),row.names=FALSE,quote=FALSE,sep="\t")
}