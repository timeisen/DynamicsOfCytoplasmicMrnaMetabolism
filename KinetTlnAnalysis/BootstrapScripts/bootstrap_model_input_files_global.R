#This script generates the single tails for model inputs
library(plyr)

for(i in 1:10){
	print("Reading background subtracted files...")
	#keep steady state in there?
	tp40m_ <- read.table(paste0("single_tag_files/background_subtracted_single_tag_files/bootstrap/tp40m_bs_minus_miR-155_v6_BOOTSTRAP_",i,".txt"))
	tp1hr_ <- read.table(paste0("single_tag_files/background_subtracted_single_tag_files/bootstrap/tp1hr_bs_minus_miR-155_v6_BOOTSTRAP_",i,".txt"))
	tp2hr_ <- read.table(paste0("single_tag_files/background_subtracted_single_tag_files/bootstrap/tp2hr_bs_minus_miR-155_v6_BOOTSTRAP_",i,".txt"))
	tp4hr_ <- read.table(paste0("single_tag_files/background_subtracted_single_tag_files/bootstrap/tp4hr_bs_minus_miR-155_v6_BOOTSTRAP_",i,".txt"))
	tp8hr_ <- read.table(paste0("single_tag_files/background_subtracted_single_tag_files/bootstrap/tp8hr_bs_minus_miR-155_v6_BOOTSTRAP_",i,".txt"))
	SS____ <- read.table(paste0("single_tag_files/background_subtracted_single_tag_files/bootstrap/tpSS_minus_norm_miR-155_BOOTSTRAP_",i,".txt"))
	
	# SS____[,253] <- SS____[,1]
	# SS____[,1] <- NULL
	
	#remove the 250
	tp40m_<-tp40m_[,-251]
	tp1hr_<-tp1hr_[,-251]
	tp2hr_<-tp2hr_[,-251]
	tp4hr_<-tp4hr_[,-251]
	tp8hr_<-tp8hr_[,-251]
	SS____<-SS____[,-251]
	
	colnames(tp40m_)<-c(paste0("tp40m_",0:249),"accession")
	colnames(tp1hr_)<-c(paste0("tp1hr_",0:249),"accession")
	colnames(tp2hr_)<-c(paste0("tp2hr_",0:249),"accession")
	colnames(tp4hr_)<-c(paste0("tp4hr_",0:249),"accession")
	colnames(tp8hr_)<-c(paste0("tp8hr_",0:249),"accession")
	colnames(SS____)<-c(paste0("SS____",0:249),"accession")
	
	tp40m_ <- tp40m_[,251:1]
	tp1hr_ <- tp1hr_[,251:1]
	tp2hr_ <- tp2hr_[,251:1]
	tp4hr_ <- tp4hr_[,251:1]
	tp8hr_ <- tp8hr_[,251:1]
	SS____ <- SS____[,251:1]
	
	m<-join_all(list(
		tp40m_,
		tp1hr_,
		tp2hr_,
		tp4hr_,
		tp8hr_,
		SS____),by="accession",type="full")
	
	m<-m[complete.cases(m),]
	m <- m[sample(nrow(m)),] #permute rows
	maxrows <- floor(nrow(m)/100)
	print("Writing permuted 100 gene files...")

	for (j in 1:maxrows){
		mi<-m[(1+(100)*(j-1)):(100*j),] #100 rows
		
		mt<-rbind(
			c(t(mi[,2:251])),
			c(t(mi[,252:501])),
			c(t(mi[,502:751])),
			c(t(mi[,752:1001])),
			c(t(mi[,1002:1251])),
			c(t(mi[,1252:1501])))
		colnames(mt)<-rep(mi$accession,each = 250)
	
		fn = paste0("model_input_files/randomized_gene_sets/bootstrap/miR-155_minus_rand_g_",j,"_BOOTSTRAP_",i,".txt")
		write.table(mt,file=fn,row.names=FALSE,quote=FALSE,sep="\t")
	}	
}
	