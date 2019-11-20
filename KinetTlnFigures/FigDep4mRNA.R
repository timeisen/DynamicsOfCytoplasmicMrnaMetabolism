library(plyr)
library(reshape2)
library(ggplot2)
source("ggplot_theme.R")

a <- list(read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/single_tag_files/background_subtracted_single_tag_files/SS_minus_norm_miR-1_v7_HYBRID20190708.txt"),read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/single_tag_files/background_subtracted_single_tag_files/SS_minus_norm_miR-155_v7_HYBRID20190708.txt"))
dat = NULL

metagene <- function(df_orig){
	df <- df_orig
	df$rs <- rowSums(df)
	m = median(df$rs)
	df[which(df$rs > m),-252] <- t(apply(df[which(df$rs > m),],1,function(x,m){x[1:251]/x[252]*m},m = m))
	return(df)

}
for(i in c(1:2)){
	ldf<-NULL
	sdf<-NULL
	miR1_SS <- a[[i]]


	# hl<-read.table("/lab/solexa_bartel/eichhorn/5EU_RPF/miR-1_miR-155_analysis/Halflife_analysis/Halflife_measurements/	miR-1_RNA_minusdox_halflife_SS.txt",head=TRUE,sep="\t")
	# hl <- hl[,c(1,7,8)]
	# hl$halflife <- log(2)/hl$beta
	# hl <- hl[,c(1,4)]
	hl<-read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/processed_files/halflifeComparisons/miR-155_minu_PAL_halflives_logspace_global_offset_SS_RNAseqScaling.txt",head=TRUE,sep="\t")
	colnames(miR1_SS) <- c(250:0,"accession")
	
	# colnames(miR1_SS) <- c("accession",0:250)
	miR1_SS <- merge(miR1_SS,hl,by="accession")
	
	miR1_SS <- miR1_SS[order(-miR1_SS$halflife),]
	bin_no <- 30
	miR1_SS$bin <- c(rep(1:bin_no,each = floor(nrow(miR1_SS)/bin_no)),rep(NA,nrow(miR1_SS)-floor(nrow(miR1_SS)/bin_no)*bin_no))
	# miR1_SS$mean <- apply(miR1_SS,1,function(x){weighted.mean(0:250,w = x[2:252])})
	# miR1_SS$sum <- rowSums(miR1_SS[,2:252])

	# miR1_bins <- data.frame(t(data.frame(lapply(1:bin_no,function(x){colSums(miR1_SS[which(miR1_SS$bin==x),2:252])/sum(miR1_SS[which(miR1_SS$bin==x),2:252])}))))
	# rownames(miR1_bins)<-NULL
	# miR1_bins$bin <- as.factor(1:bin_no)
	# colnames(miR1_bins)<-c(0:250,"bin")
	
	# miR1_melt <- melt(miR1_bins,id="bin",value.name="abundance",variable.name="tail_length")
	# miR1_melt$tail_length <- as.numeric(as.character(miR1_melt$tail_length))
	
	# p1 <- ggplot(miR1_melt,aes(x=tail_length,y=abundance))+geom_line()+facet_wrap(~bin)+theme_bw()+geom_smooth(span=0.1)

	#for hl_bin
	long = 9
	short = .5
	miR1_SS_thresh <- miR1_SS[which(miR1_SS$halflife > long | miR1_SS$halflife < short),]
	miR1_SS_thresh[which(miR1_SS_thresh$halflife > long),]$bin = "Long"
	miR1_SS_thresh[which(miR1_SS_thresh$halflife < short),]$bin = "Short"
	
	# ldf = metagene(miR1_SS_thresh[which(miR1_SS_thresh$bin=="Long"),-c(1,253:254)])
	# sdf = metagene(miR1_SS_thresh[which(miR1_SS_thresh$bin=="Short"),-c(1,253:254)])

	ldf = miR1_SS_thresh[which(miR1_SS_thresh$bin=="Long"),-c(1,253:254)]
	sdf = miR1_SS_thresh[which(miR1_SS_thresh$bin=="Short"),-c(1,253:254)]

	ldf$bin = "Long"
	sdf$bin = "Short"

	miR1_SS_thresh = rbind(ldf,sdf)
	print(nrow(ldf))
	print(nrow(sdf))

	miR1_SS_thresh_bin <- data.frame(t(data.frame(lapply(c("Long","Short"),function(x){colSums(miR1_SS_thresh[which(miR1_SS_thresh$bin==x),1:250])/sum(	miR1_SS_thresh[which(miR1_SS_thresh$bin==x),1:250])}))))
	rownames(miR1_SS_thresh_bin)<-NULL
	miR1_SS_thresh_bin$bin <- c("Long","Short")
	colnames(miR1_SS_thresh_bin)<-c(0:249,"bin")
	
	
	miR1_SS_thresh_bin_melt <- melt(miR1_SS_thresh_bin,id="bin",value.name="abundance",variable.name="tail_length")
	miR1_SS_thresh_bin_melt$tail_length <- as.numeric(as.character(miR1_SS_thresh_bin_melt$tail_length))
	miR1_SS_thresh_bin_melt$sample <- c("miR-1","miR-155")[i]
	dat <- rbind(dat,miR1_SS_thresh_bin_melt)
}

#plotting

p2 <- ggplot(dat,aes(x=tail_length,y=abundance))+geom_line()+facet_wrap(~sample+bin)+theme_tim()+geom_smooth(data=dat,aes(x=tail_length,y=abundance),span=0.1,method="loess")+scale_y_continuous(limits=c(0,0.01),name=NULL,expand=c(0,0))+scale_x_continuous(expand=c(0,0))
p3a <- ggplot(dat[which(dat$sample == "miR-155" & dat$bin == "Long"),],aes(x=tail_length,y=abundance)) + 
	  geom_line() +  
	  theme_tim() + 
	  #geom_smooth(data=dat,aes(x=tail_length,y=abundance),span=0.1,method="loess") + 
	  scale_y_continuous(limits=c(0,0.01),name=NULL,expand=c(0,0))+scale_x_continuous(expand=c(0,0),limits=c(0,250))

p3b <- ggplot(dat[which(dat$sample == "miR-155" & dat$bin == "Short"),],aes(x=tail_length,y=abundance)) + 
	  geom_line() + 
	  theme_tim() + 
	  #geom_smooth(data=dat,aes(x=tail_length,y=abundance),span=0.1,method="loess") + 
	  scale_y_continuous(limits=c(0,0.01),name=NULL,expand=c(0,0))+scale_x_continuous(expand=c(0,0),limits=c(0,250))

dat155 <- dat[which(dat$sample == "miR-155"),]
dat155$bin <- factor(dat155$bin, levels = c("Short","Long")) #plot short first

p3 <- ggplot(dat155,aes(x = tail_length,y = abundance, color = bin)) + 
	  geom_step() +  
	  theme_tim() + 
	  #geom_smooth(data=dat,aes(x=tail_length,y=abundance),span=0.1,method="loess") + 
	  scale_y_continuous(limits=c(0,0.01),name=NULL,expand=c(0,0)) + 
	  scale_x_continuous(expand=c(0,0),limits=c(0,250)) +
	  scale_color_manual(values = c("blue","red"))

p4 <- ggplot(dat155,aes(x = tail_length,y = abundance, color = bin)) + 
	  geom_step() +  
	  theme_tim() + 
	  #geom_smooth(data=dat,aes(x=tail_length,y=abundance),span=0.1,method="loess") + 
	  scale_y_continuous(limits=c(0,0.01),name=NULL,expand=c(0,0)) + 
	  scale_x_continuous(expand=c(0,0),limits=c(0,250)) +
	  scale_color_manual(values = c("blue","red"))

ggsave(plot=p3,file="figures/version_8/FigS8A.pdf",width=2,height=2)


# break 
# #periodicity analysis
# test <- dat[which(dat$bin=="Long" & dat$sample=="miR-155"),]
# test <- test[-c(1:5,251),]

# period <- function(test){
# 	lf = predict(loess(test$abundance~test$tail_length))
# 	p2 = periodogram(test$abundance/lf)
# 	dd2 = data.frame(freq=p2$freq, spec=p2$spec)
# 	dd2$time <- 1/dd2$freq
# 	dd2 = dd2[order(-dd2$spec),]
# 	return(dd2[1:10,3])}

# print(period(test))



