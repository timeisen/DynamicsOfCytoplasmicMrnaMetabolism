library(ggplot2)
library(plyr)
library(reshape2)
options("scipen"=5)
source("/lab/solexa_bartel/teisen/RNAseq/Scripts/general/ggplot_theme.R")

stds <- read.table("/lab/solexa_bartel/teisen/RNAseq/Annotation_files/TAIL_seq/list_of_standards_SWE_notation.txt")$V1
hl <- read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/processed_files/halflifeComparisons/miR-155_minu_PAL_halflives_logspace_global_offset_SS.txt",head=TRUE)
mSWE <- read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/other_analyses/directLig/SE_column_3T3_cDNA_ALL_READS_18_08_17/AATCCG-s_8_1_sequence_singletails_pairedEndMapping_reformat_50tagCutoff_ALL.txt") #SWE direct 3p Lig dataset
mSIM <- read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/other_analyses/all_genes_simulated_V75.txt")
mSIM[2:251] <- mSIM[251:2]
dead <- read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/rate_constant_measurements/miR-155_minus_samples_UNLINKV75_run6_final_global_param.txt")
dead <- dead[,c(1,4,5)]
colnames(dead)<-c("accession","dead_rate","dcp_rate")

colnames(mSWE)[1] <- "accession"
mPAL <- read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/single_tag_files/background_subtracted_single_tag_files/SS_minus_norm_miR-155_v6.txt")
# mPAL <- read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/single_tag_files/background_subtracted_single_tag_files/SS_minus_miR-155_v5.txt")

colnames(mPAL)[252] <- "accession"
# colnames(mPAL)[1] <- "accession"
colnames(mSIM)[1] <- "accession"
mSWE <- mSWE[which(!mSWE$accession %in% stds),]
mPAL <- mPAL[which(!mPAL$accession %in% stds),]
mSIM <- mSIM[which(!mSIM$accession %in% stds),]

mSWEmod <- mSWE

mSWEmod$rs = rowSums(mSWEmod[,-1])
mSWEaccession = as.character(mSWEmod[,1])
mSWEmod = mSWEmod[,2:9]/mSWEmod$rs
mSWEmod$accession = mSWEaccession
mPAL$accession = as.character(mPAL$accession)
mPALSWE <- merge(mPAL,mSWEmod,by="accession")
# break #these lines need to be fixed
mPALSWE[,2:9] = rowSums(mPALSWE[,2:252])*mPALSWE[,253:260]
mPALSWE[,253:260] <- NULL

mPALSWE <- mPALSWE[which(!mPALSWE$accession %in% stds),]
mPALSWE_write <- mPALSWE[,c(2:ncol(mPALSWE),1)]
# write.table(mPALSWE_write,file="/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/single_tag_files/background_subtracted_single_tag_files/SS_minus_norm_miR-155_v6_HYBRID.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
csPAL <- colSums(mPAL[mPAL$accession %in% mSWE$accession,][,1:248])/sum(mPAL[mPAL$accession %in% mSWE$accession,-252])
csSWE <- colSums(mSWE[mSWE$accession %in% mPAL$accession,][,2:249])/sum(mSWE[mSWE$accession %in% mPAL$accession,-1])
mSWE <- mSWE[which(!mSWE$accession %in% stds),]
mPAL <- mPAL[which(!mPAL$accession %in% stds),]
mSIM <- mSIM[which(!mSIM$accession %in% stds),]

mPALavg <- data.frame(
	"accession" = mPAL$accession, 
	"mPALavg" = apply(mPAL[,1:8],1,mean)/apply(mPAL[,51:151],1,mean))
mSWEavg <- data.frame(
	"accession" = mSWE$accession, 
	"mSWEavg" = apply(mSWE[,2:9],1,mean)/apply(mSWE[,52:152],1,mean))
mPALSWEavg <- data.frame(
	"accession" = mPALSWE$accession, 
	"mPALSWEavg" = apply(mPALSWE[,2:9],1,mean)/apply(mPALSWE[,52:152],1,mean))
mavg <- join_all(list(mPALavg,mSWEavg,mPALSWEavg,dead),type="full",by="accession")
mavg <- mavg[complete.cases(mavg),]
mavgMelt <- melt(mavg,id=c("dead_rate","dcp_rate","accession"))
annot<-read.table("/lab/solexa_bartel/teisen/RNAseq/Annotation_files/TAIL_seq/symbol_to_accession.txt",head=TRUE)
mavgMelt<-merge(mavgMelt,annot,by="accession")
pAvg <- ggplot(mavgMelt,aes(x=value,y=dead_rate,labels=symbol))+geom_point(size=0.5)+theme_bw()+facet_wrap(~variable)+geom_abline()+scale_x_continuous(trans="log10",name="mean of < 8nt / mean of 50-250")+scale_y_continuous(trans="log10")
# mSWE$V2 = 0
# mSWE$V2 = 0
# mSWE$V3 = 0
# mSWE$V4 = 0
# mSWE$V5 = 0
# mSWE$V6 = 0
# mSWE$V7 = 0
# mSWE$V8 = 0

# mSWE <- head(mSWE)
mSWE_frac <- data.frame(sapply(2:252,function(x){
	if(x==2){return(mSWE[,2]/rowSums(mSWE[,-1]))}
	rowSums(mSWE[,2:x])/rowSums(mSWE[,-1])}))
mSWE_frac$accession <- mSWE$accession
mSWE_frac <- merge(mSWE_frac,hl[,c(1,5)],by="accession")


mPALSWE_frac <- data.frame(sapply(2:252,function(x){
	if(x==2){return(mPALSWE[,2]/rowSums(mPALSWE[,-1]))}
	rowSums(mPALSWE[,2:x])/rowSums(mPALSWE[,-1])}))
mPALSWE_frac$accession <- mPALSWE$accession
mPALSWE_frac <- merge(mPALSWE_frac,hl[,c(1,5)],by="accession")

mPALSWE_rev <- mPALSWE[,ncol(mPALSWE):1]
mPALSWE_frac_rev <- data.frame(sapply(1:251,function(x){
	if(x==1){return(mPALSWE_rev[,1]/rowSums(mPALSWE_rev[,-252]))}
	rowSums(mPALSWE_rev[,1:x])/rowSums(mPALSWE_rev[,-252])}))
mPALSWE_frac_rev$accession <- mPALSWE_rev$accession
mPALSWE_frac_rev <- merge(mPALSWE_frac_rev,hl[,c(1,5)],by="accession")

# mPALSWErev <- mPALSWE#[,ncol(mPALSWE):1]
# mPALSWErev_frac <- data.frame(sapply(2:242,function(x){
# 	rowSums(mPALSWErev[,x:(x+10)])/rowSums(mPALSWErev[,-1])}))
# mPALSWErev_frac$accession <- mPALSWErev$accession
# mPALSWErev_frac <- merge(mPALSWErev_frac,hl[,c(1,5)],by="accession")

mPAL_frac <- data.frame(sapply(1:251,function(x){
	if(x==1){return(mPAL[,1]/rowSums(mPAL[,-252]))}
	rowSums(mPAL[,1:x])/rowSums(mPAL[,-252])}))
mPAL_frac$accession <- mPAL$accession
mPAL_frac <- merge(mPAL_frac,hl[,c(1,5)],by="accession")

# mSIM <- head(mSIM)
mSIM_frac <- data.frame(sapply(2:251,function(x){
	if(x==2){return(mSIM[,2]/rowSums(mSIM[,-1]))}
	rowSums(mSIM[,2:x])/rowSums(mSIM[,-1])}))
mSIM_frac$accession <- mSIM$accession
mSIM_frac <- merge(mSIM_frac,hl[,c(1,5)],by="accession")
# break
# corData <- data.frame(
# 	tl = 1:251,
# 	mPALSWE_cor = as.numeric(tail(cor(log(mPALSWE_frac[,-1]+.001),method="pearson",use="complete.obs"),1))[-252],
# 	# mPALSWE_cor_rev = as.numeric(tail(cor(mPALSWErev_frac[,-1],method="spearman"),1))[-252]) 
# 	mSWE_cor = as.numeric(tail(cor(log(mSWE_frac[,-1]+.001),method="pearson",use="complete.obs"),1))[-252],
# 	mPAL_cor = as.numeric(tail(cor(log(mPAL_frac[,-1]+.001),method="pearson",use="complete.obs"),1))[-252])
# 	mSIM_cor = c(as.numeric(tail(cor(mSIM_frac[,-1],method="spearman"),1)[-251]),NA))

corData <- data.frame(
	tl = 1:251,
	mPALSWE_cor = as.numeric(tail(cor(mPALSWE_frac[,-1],method="spearman",use="complete.obs"),1))[-252],
	mSWE_cor = as.numeric(tail(cor(mSWE_frac[,-1],method="spearman",use="complete.obs"),1))[-252],
	mPALSWE_cor_rev = as.numeric(tail(cor(mPALSWE_frac_rev[,-1],method="spearman"),1))[-252],
	mPAL_cor = as.numeric(tail(cor(mPAL_frac[,-1],method="spearman",use="complete.obs"),1))[-252])

# corDataBin <- data.frame(
# 	tl = 1:242,
# 	mPALSWE_cor_rev = as.numeric(tail(cor(mPALSWErev_frac[,-1],method="spearman"),1))[-252])

corDataMelt <- melt(corData,id="tl")
# corDataBinMelt <- melt(corDataBin,id="tl")

p1 <- ggplot(corDataMelt[which(corDataMelt$variable=="mPALSWE_cor"),],aes(x=tl,y=value,color=variable))+geom_step()+theme_tim()+scale_y_continuous(name=NULL,expand=c(0,0),limits=c(-1,1))+scale_x_continuous(expand=c(0,0),limits=c(0,250))
p1a <- ggplot(corDataMelt[which(corDataMelt$variable=="mPALSWE_cor"),],aes(x=tl,y=value,color=variable))+geom_step()+theme_tim()+scale_y_continuous(name=NULL,expand=c(0,0),limits=c(-1,1))+scale_x_continuous(expand=c(0,0),limits=c(0,100))
p1rev <- ggplot(corDataMelt[which(corDataMelt$variable=="mPALSWE_cor_rev"),],aes(x=tl,y=value,color=variable))+geom_step()+theme_tim()+scale_y_continuous(name=NULL,expand=c(0,0),limits=c(-1,1))+scale_x_continuous(expand=c(0,0),limits=c(0,150),breaks=c(0:6)*25,labels=(10:4)*25)

#Now the scatter plots ....
print(nrow(mPALSWE_frac))
print(nrow(mPAL_frac))
print(nrow(mSWE_frac))
print(cor(mPALSWE_frac$halflife_t,mPALSWE_frac$X19,method="spearman",use="complete.obs"))
print(cor(mPAL_frac$halflife_t,mPAL_frac$X19,method="spearman"))
print(cor(mSWE_frac$halflife_t,mSWE_frac$X19,method="spearman"))
print(cor(mPALSWE_frac_rev$halflife_t,mPALSWE_frac_rev$X75,method="spearman"))

weightedCV = function(x){
	mn = weighted.mean(0:250,w = x)
	sdev = (sum(x * (0:250 - mn)^2)/sum(x))^0.5
	#returning only the sdev here
	return(sdev)}

mPALSWEdisp <- mPALSWE
mPALSWEdisp$cv <- apply(mPALSWEdisp[,-1],1,weightedCV)
mPALSWEdisp <- merge(mPALSWEdisp,hl[,c(1,5)],by="accession")
print(cor(mPALSWEdisp$halflife_t,mPALSWEdisp$cv,method="spearman"))

p5cv <- ggplot(mPALSWEdisp,aes(x=halflife_t,y=cv))+geom_point(alpha=0.5,size=0.5)+theme_tim()+scale_x_continuous(trans="log10",breaks=log_ticks(0.01,100)[[2]],labels=log_ticks(0.01,100)[[1]],expand=c(0,0))+scale_y_continuous(name=NULL,expand=c(0,0),limits=c(0,100))
ggsave(plot=p5cv,width=2.0,height=2.0,useDingbats=FALSE,file="/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/figures/other/FigS11Bpanelcv.pdf")

# p2 <- ggplot(mPAL_frac,aes(x=halflife_t,y=X19))+geom_point(alpha=0.5,size=0.5)+theme_tim()+scale_x_continuous(trans="log10",name=NULL,limits=c(0.01,100),breaks=log_ticks(0.01,100)[[2]],labels=log_ticks(0.01,100)[[1]],expand=c(0,0))+scale_y_continuous(name=NULL,breaks=log_ticks(0.0001,1)[[2]],labels=log_ticks(0.0001,1)[[1]],expand=c(0,0),trans="log10",limits=c(0.0001,1))
# p3 <- ggplot(mSWE_frac,aes(x=halflife_t,y=X19))+geom_point(alpha=0.5,size=0.5)+theme_tim()+scale_x_continuous(trans="log10",breaks=log_ticks(0.01,100)[[2]],labels=log_ticks(0.01,100)[[1]],expand=c(0,0))+scale_y_continuous(name=NULL,breaks=log_ticks(0.0001,1)[[2]],labels=log_ticks(0.0001,1)[[1]],expand=c(0,0),trans="log10",limits=c(0.0001,1))
mPALSWE_frac[mPALSWE_frac$halflife_t>100,]$halflife_t = 100
# mPALSWE_frac[mPALSWE_frac$X19==0,]$X19 = 0.001
p4 <- ggplot(mPALSWE_frac,aes(x=halflife_t,y=X19))+geom_point(alpha=0.5,size=0.5)+theme_tim()+scale_x_continuous(trans="log10",breaks=log_ticks(0.01,100)[[2]],labels=log_ticks(0.01,100)[[1]],expand=c(0,0))+scale_y_continuous(name=NULL,breaks=log_ticks(0.0001,1)[[2]],labels=log_ticks(0.0001,1)[[1]],expand=c(0,0),trans="log10",limits=c(0.0001,1))


# this is for the large fraction correlation
# p5 <- ggplot(mPALSWE_frac,aes(x=halflife_t,y=X172))+geom_point(alpha=0.5,size=0.5)+theme_tim()+scale_x_continuous(trans="log10",breaks=log_ticks(0.01,100)[[2]],labels=log_ticks(0.01,100)[[1]],expand=c(0,0))+scale_y_continuous(name=NULL,breaks=log_ticks(0.0001,1)[[2]],labels=log_ticks(0.0001,1)[[1]],expand=c(0,0),trans="log10",limits=c(0.0001,1))

p5b <- ggplot(mPALSWE_frac_rev,aes(x=halflife_t,y=X75))+geom_point(alpha=0.5,size=0.5)+theme_tim()+scale_x_continuous(trans="log10",breaks=log_ticks(0.01,100)[[2]],labels=log_ticks(0.01,100)[[1]],expand=c(0,0))+scale_y_continuous(expand=c(0,0),limits=c(0.01,1),trans="log10",breaks=log_ticks(0.001,1)[[2]],labels=log_ticks(0.001,1)[[1]],name=NULL)
# break
# p4 <- ggplot(mSIM_frac,aes(x=halflife_t,y=X19))+geom_point(alpha=0.5,size=0.5)+theme_tim()+scale_x_continuous(trans="log10",breaks=log_ticks(0.01,100)[[2]],labels=log_ticks(0.01,100)[[1]],expand=c(0,0))+scale_y_continuous(name=NULL,breaks=log_ticks(0.0001,1)[[2]],labels=log_ticks(0.0001,1)[[1]],expand=c(0,0),trans="log10",limits=c(0.0001,1))
# print(cor(mSIM_frac$halflife_t,mSIM_frac$X19))

# ggsave(plot=p4,width=2,height=2,useDingbats=FALSE,file="figures/other/buildupIndividlDcpModel.pdf")

ggsave(plot=p5b,width=2.0,height=2.0,useDingbats=FALSE,file="/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/figures/other/FigS11BpanelLongTails.pdf")
ggsave(plot=p1rev,width=2.0,height=2.0,useDingbats=FALSE,file="/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/figures/other/FigS11CpanelLongTails.pdf")

# ggsave(plot=p1,width=1.5,height=1.5,useDingbats=FALSE,file="/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/figures/other/FigS11A.pdf")
# ggsave(plot=p2,width=1.5,height=1.5,useDingbats=FALSE,file="/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/figures/other/FigS11Bpanel1.pdf")
# ggsave(plot=p3,width=1.5,height=1.5,useDingbats=FALSE,file="/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/figures/other/FigS11Bpanel2.pdf")
ggsave(plot=p4,width=2.0,height=2.0,useDingbats=FALSE,file="/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/figures/other/FigS11Bpanel3.pdf")
# ggsave(plot=p5,width=1.5,height=1.5,useDingbats=FALSE,file="/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/figures/other/FigS11BpanelX.pdf")
# ggsave(plot=p5b,width=2.0,height=2.0,useDingbats=FALSE,file="/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/figures/other/FigS11BpanelXb.pdf")
ggsave(plot=p1a,width=2.0,height=2.0,useDingbats=FALSE,file="/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/figures/other/FigS11Bpanel11A1.pdf")

# ggsave(plot=p3,width=1.5,height=1.5,useDingbats=FALSE,file="/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/figures/other/FigS11Bpanel2.pdf")
# ggsave(plot=p4,width=1.5,height=1.5,useDingbats=FALSE,file="/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/figures/other/FigS11Bpanel3.pdf")