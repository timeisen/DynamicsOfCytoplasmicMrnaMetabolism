#FigSDEP, model supplementary figure, no longer used.

library(tidyverse)
# library(ggplot2)
options("scipen"=-3)
source("/lab/solexa_bartel/teisen/RNAseq/Scripts/general/ggplot_theme.R")
#plotting parameters here
alpha = 0.1
size  = 0.5


mS <- read.table("rate_constant_measurements/miR-155_minus_samples_UNLINKV75_GlobalSt.txt")
mD <- read.table("rate_constant_measurements/miR-155_minus_samples_UNLINKV75_TwoDea.txt")
m <- read.table("rate_constant_measurements/miR-155_minus_samples_UNLINKV75_run6_final_global_param.txt")
colnames(m)<-c("accession","st","am","km","bm","rm")
colnames(mS)<-c("accession","aS","kS","bS","rS")
colnames(mD)<-c("accession","stD","aD","k1D","k2D","bD","rD")
colnames(m)<-c("accession","stm","am","km","bm","rm")
library(plyr)
mALL <- join_all(list(m,mS,mD),type="full",by="accession")

data<-read.table("model_input_files/miR-155_minus_sample_mean_tails_50tags_with_expr_background_subtracted_v6_st.txt",head=TRUE)

vars<-c(      #miR-155 minus
  2.58757e-16,
  3.053586e-16,
  2.92142e-15,
  1.16722e-14,
  2.772866e-14,
  2.394499e-13)
rss <- data.frame(
  "accession" = data[,1],
  "rss" = apply(data[,-1],1,function(x){
    sqs<-(x-mean(x))^2
    sqs[1:250]<-sqs[1:250]/vars[1]
    sqs[251:500]<-sqs[251:500]/vars[2]
    sqs[501:750]<-sqs[501:750]/vars[3]
    sqs[751:1000]<-sqs[751:1000]/vars[4]
    sqs[1001:1250]<-sqs[1001:1250]/vars[5]
    sqs[1251:1500]<-sqs[1251:1500]/vars[6]
    rss = sum(sqs*1E6)
    })
  )

twoDea <- merge(m[,c(1,6)],mD[,c(1,7)],by="accession")
st <- merge(m[,c(1,6)],mS[,c(1,5)],by="accession")
mThree <- merge(st,twoDea)
mThree<-merge(mThree,rss,by="accession")
mThree$rm<-1-mThree$rm/mThree$rss
mThree$rS<-1-mThree$rS/mThree$rss
mThree$rD<-1-mThree$rD/mThree$rss
mThree$rss <- NULL

twoDeaMelt <- gather(twoDea,key="Model",value="Resid",-accession)
stMelt <- gather(st,key="Model",value="Resid",-accession)
threeMelt <-  gather(mThree,key="Model",value="Rsq",-accession)
mBar <- data.frame(
	data = c("rm","rS","rD"),
	values = apply(mThree[,-1],2,median))

p6all <- ggplot(threeMelt,aes(x=Rsq,color=Model))+stat_ecdf()+scale_x_continuous(trans="log10",labels=fancy_scientific,breaks=log_ticks(1E-6,1)[[2]],expand=c(0,0))+theme_tim()+scale_y_continuous(name=NULL,expand=c(0,0),limits=c(0,1))+theme(axis.text.x=element_text(angle=45, hjust=1))

p7all <- ggplot(threeMelt,aes(x=Rsq,color=Model))+stat_ecdf()+scale_x_continuous(limits=c(0,1),expand=c(0,0))+coord_cartesian(xlim=c(0.5,1))+theme_tim()+scale_y_continuous(name=NULL,expand=c(0,0),limits=c(0,1))

p7bar <- ggplot(mBar,aes(x=data,y=values))+geom_bar(stat="identity")+theme_tim()+scale_y_continuous(name=NULL,expand=c(0,0))
p7box <- ggplot(threeMelt,aes(x=Model,y=Rsq))+geom_boxplot()+theme_tim()+scale_y_continuous(name=NULL,expand=c(0,0))

ggsave(plot=p7all,file="figures/version_6/FigS5A.pdf",width=2,height=2,useDingbats=FALSE)


