#Panels for figure 2C and D
library(ggplot2)
library(plyr)
library(reshape2)
library(spatstat)
# This script is an updated version of the original 2C V4 but with a background correction term. 
# 2806 genes
source("/lab/solexa_bartel/teisen/RNAseq/Scripts/general/ggplot_theme.R")
options("warn"=1)
print("Read files in ...")

#Add SS
stds<-read.table("/lab/solexa_bartel/teisen/RNAseq/Annotation_files/TAIL_seq/list_of_standards_SWE_notation.txt")

tp40m_bs<-read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/single_tag_files/background_subtracted_single_tag_files/tp40m_bs_minus_miR-155_v7.txt",sep="\t")
tp1hr_bs<-read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/single_tag_files/background_subtracted_single_tag_files/tp1hr_bs_minus_miR-155_v7.txt",sep="\t")
tp2hr_bs<-read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/single_tag_files/background_subtracted_single_tag_files/tp2hr_bs_minus_miR-155_v7.txt",sep="\t")
tp4hr_bs<-read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/single_tag_files/background_subtracted_single_tag_files/tp4hr_bs_minus_miR-155_v7.txt",sep="\t")
tp8hr_bs<-read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/single_tag_files/background_subtracted_single_tag_files/tp8hr_bs_minus_miR-155_v7.txt",sep="\t")
tpSS_bs<-read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/single_tag_files/background_subtracted_single_tag_files/SS_minus_norm_miR-155_v7_SsScaled.txt",sep="\t")


tp40m_bs_norm <- sum(tp40m_bs[,-c(252)])
tp1hr_bs_norm <- sum(tp1hr_bs[,-c(252)])
tp2hr_bs_norm <- sum(tp2hr_bs[,-c(252)])
tp4hr_bs_norm <- sum(tp4hr_bs[,-c(252)])
tp8hr_bs_norm <- sum(tp8hr_bs[,-c(252)])
 tpSS_bs_norm <- sum( tpSS_bs[,-c(252)])

colnames(tp40m_bs)[252] <- "accession"
colnames(tp1hr_bs)[252] <- "accession"
colnames(tp2hr_bs)[252] <- "accession"
colnames(tp4hr_bs)[252] <- "accession"
colnames(tp8hr_bs)[252] <- "accession"
colnames(tpSS_bs)[252] <- "accession"

############## Added 2019 05 16 ##############
# tp40m_bs_t <- tp40m_bs[,-252]
# tp1hr_bs_t <- tp1hr_bs[,-252]
# tp2hr_bs_t <- tp2hr_bs[,-252]
# tp4hr_bs_t <- tp4hr_bs[,-252]
# tp8hr_bs_t <- tp8hr_bs[,-252]
#  tpSS_bs_t <-  tpSS_bs[,-252]

# tp40m_bs_t[tp40m_bs_t < 0] = 0 
# tp1hr_bs_t[tp1hr_bs_t < 0] = 0 
# tp2hr_bs_t[tp2hr_bs_t < 0] = 0 
# tp4hr_bs_t[tp4hr_bs_t < 0] = 0 
# tp8hr_bs_t[tp8hr_bs_t < 0] = 0 
#  tpSS_bs_t[ tpSS_bs_t < 0] = 0 

# tp40m_bs[,-252] <- tp40m_bs_t
# tp1hr_bs[,-252] <- tp1hr_bs_t
# tp2hr_bs[,-252] <- tp2hr_bs_t
# tp4hr_bs[,-252] <- tp4hr_bs_t
# tp8hr_bs[,-252] <- tp8hr_bs_t
#  tpSS_bs[,-252] <-  tpSS_bs_t
############## Added 2019 05 16 ##############


tp40m_bs_mean<-data.frame(apply(tp40m_bs[,-252],1,function(x) weighted.mean(0:250,w=x)))
tp1hr_bs_mean<-data.frame(apply(tp1hr_bs[,-252],1,function(x) weighted.mean(0:250,w=x)))
tp2hr_bs_mean<-data.frame(apply(tp2hr_bs[,-252],1,function(x) weighted.mean(0:250,w=x)))
tp4hr_bs_mean<-data.frame(apply(tp4hr_bs[,-252],1,function(x) weighted.mean(0:250,w=x)))
tp8hr_bs_mean<-data.frame(apply(tp8hr_bs[,-252],1,function(x) weighted.mean(0:250,w=x)))
tpSS_mean<-data.frame(apply(tpSS_bs[,-252],1,function(x) weighted.mean(0:250,w=x)))


tp40m_bs_mean$accession<-tp40m_bs[,252]
tp1hr_bs_mean$accession<-tp1hr_bs[,252]
tp2hr_bs_mean$accession<-tp2hr_bs[,252]
tp4hr_bs_mean$accession<-tp4hr_bs[,252]
tp8hr_bs_mean$accession<-tp8hr_bs[,252]
tpSS_mean$accession<-tpSS_bs[,252]


bins<-function(df,nt_bins){
	df<-df[order(df[,2]),]
	x<-cut(df[,2],seq(0,250,by=nt_bins),seq(1,250/nt_bins,by=1),include.lowest = TRUE)
	df$x<-x
	x<-unique(x)
	df2<-data.frame(sapply(1:length(x),function(n){return(sum(df[which(df$x==x[n]),][,1]))}),x)
	return(df2)}
bw=2

tp40m_bs <- colSums(tp40m_bs[,-c(251,252)])
tp1hr_bs <- colSums(tp1hr_bs[,-c(251,252)])
tp2hr_bs <- colSums(tp2hr_bs[,-c(251,252)])
tp4hr_bs <- colSums(tp4hr_bs[,-c(251,252)])
tp8hr_bs <- colSums(tp8hr_bs[,-c(251,252)])
tpSS_bs  <- colSums(tpSS_bs[,-c(251,252)])

tp40m_bs <- tp40m_bs/sum(tp40m_bs)
tp1hr_bs <- tp1hr_bs/sum(tp1hr_bs)
tp2hr_bs <- tp2hr_bs/sum(tp2hr_bs)
tp4hr_bs <- tp4hr_bs/sum(tp4hr_bs)
tp8hr_bs <- tp8hr_bs/sum(tp8hr_bs)
tpSS_bs <- tpSS_bs/sum(tpSS_bs)

tp40m_bs<-data.frame(tp40m_bs)
tp1hr_bs<-data.frame(tp1hr_bs)
tp2hr_bs<-data.frame(tp2hr_bs)
tp4hr_bs<-data.frame(tp4hr_bs)
tp8hr_bs<-data.frame(tp8hr_bs)
tpSS_bs<-data.frame(tpSS_bs)


print("Calculating bins ...")
tp40m_bs$tail_length<-0:249
tp1hr_bs$tail_length<-0:249
tp2hr_bs$tail_length<-0:249
tp4hr_bs$tail_length<-0:249
tp8hr_bs$tail_length<-0:249
tpSS_bs$tail_length<-0:249


tp40m_bs_bw<-bins(tp40m_bs,bw)
tp1hr_bs_bw<-bins(tp1hr_bs,bw)
tp2hr_bs_bw<-bins(tp2hr_bs,bw)
tp4hr_bs_bw<-bins(tp4hr_bs,bw)
tp8hr_bs_bw<-bins(tp8hr_bs,bw)
tpSS_bs_bw<-bins(tpSS_bs,bw)

colnames(tp40m_bs_bw)<-"40 m"
colnames(tp1hr_bs_bw)<-"1 h"
colnames(tp2hr_bs_bw)<-"2 h"
colnames(tp4hr_bs_bw)<-"4 h"
colnames(tp8hr_bs_bw)<-"8 h"
colnames(tpSS_bs_bw)<-"SS"


tp40m_bs_bw[,2]<-NULL
tp1hr_bs_bw[,2]<-NULL
tp2hr_bs_bw[,2]<-NULL
tp4hr_bs_bw[,2]<-NULL
tp8hr_bs_bw[,2]<-NULL
tpSS_bs_bw[,2]<-NULL


tails<-cbind(
	tp40m_bs_bw,
	tp1hr_bs_bw,
	tp2hr_bs_bw,
	tp4hr_bs_bw,
	tp8hr_bs_bw,
	tpSS_bs_bw)

tails$tail_length<-1:(250/bw)
cols<-c(
"black",
"#8856a7",
"#2b8cbe",
"#31a354",
"#fec44f",
"#f03b20")

cols<-rev(cols)


tails_melt<-melt(tails,id="tail_length")

tails_melt$norm<-0
# expr_data <- read.table("processed_files/rnaseqDatasets/miR-155_minu_polyA_30minHLnorm.txt",head = TRUE,sep = "\t")
# expr_norm <- colSums(expr_data[,-1])
print(tp40m_bs_norm)
print(tp1hr_bs_norm)
print(tp2hr_bs_norm)
print(tp4hr_bs_norm)
print(tp8hr_bs_norm)
print(tpSS_bs_norm)

print(tpSS_bs_norm/tp8hr_bs_norm)
print(tp8hr_bs_norm/tp4hr_bs_norm)
print(tp4hr_bs_norm/tp2hr_bs_norm)
print(tp2hr_bs_norm/tp1hr_bs_norm)
print(tp1hr_bs_norm/tp40m_bs_norm)



tails_melt[which(tails_melt$variable == "40 m"),]$norm <-  1 / tp40m_bs_norm#1
tails_melt[which(tails_melt$variable == "1 h"),]$norm  <-  1 / tp1hr_bs_norm#1.74449207
tails_melt[which(tails_melt$variable == "2 h"),]$norm  <-  1 / tp2hr_bs_norm#0.193279838
tails_melt[which(tails_melt$variable == "4 h"),]$norm  <-  1 / tp4hr_bs_norm#0.12829349
tails_melt[which(tails_melt$variable == "8 h"),]$norm  <-  1 / tp8hr_bs_norm#0.079662677
tails_melt[which(tails_melt$variable == "SS"),]$norm   <-  1 / tpSS_bs_norm#0.048435329
 

#tails<-rbind(tp40min,tp1hr)
tails_melt$variable<-factor(tails_melt$variable,levels=c("40 m","1 h","2 h","4 h","8 h","SS"))

p1<-ggplot(tails_melt,aes(x=tail_length*bw,y=value,color=variable)) +
geom_step(size=0.5) +
scale_y_continuous(name=NULL,limits=c(0,0.03),expand=c(0,0)) +
theme_tim() +
scale_x_continuous(name=NULL,limits=c(8,250), breaks = c(8,(1:5) * 50), labels = c(8,(1:5) * 50)) + 
scale_colour_manual(values=cols) +
theme(axis.line.x = element_blank()) + 
geom_segment(aes(x = 8, xend = 250, y = 0, yend = 0))


# p1<-ggplot(tails,aes(x=tail_length*5,y=count/norm1,color=sample))+
# geom_step(size=0.5)+scale_y_continuous(name=NULL,limits=c(0,0.06),expand=c(0,0))+theme_tim()+
# scale_x_continuous(name=NULL,expand=c(0,3),limits=c(0,250))+scale_colour_manual(values=cols)
tails_melt$abund <- tails_melt$value/tails_melt$norm
tails_melt$abund <- tails_melt$abund/sum(tails_melt[which(tails_melt$variable == "SS"),]$abund)

p2<-ggplot(tails_melt,aes(x=tail_length*bw,y=abund,color=variable)) + 
geom_step(size=0.5) + 
scale_y_continuous(name=NULL,expand=c(0,0)) + 
theme_tim() +
scale_x_continuous(name=NULL,limits=c(8,250), breaks = c(8,(1:5) * 50), labels = c(8,(1:5) * 50)) + 
scale_colour_manual(values=cols) +
theme(axis.line.x = element_blank()) + 
geom_segment(aes(x = 8, xend = 250, y = 0, yend = 0))


# p1<-ggplot(tails,aes(x=tail_length,y=..count../sum(..count..),weights=counts,colour=sample))+
# stat_bin(geom="step",bins=50,position='identity',size=.5)+theme_classic()+
# scale_y_continuous(name=NULL)+
# scale_x_continuous(name=NULL)
# #theme(legend.position="none")

# p2<-ggplot(tails,aes(x=tail_length,y=..density..,colour=sample))+
# stat_bin(geom="step",bins=50,position='identity',size=.5)+theme_classic()+
# scale_y_continuous(name=NULL)+
# scale_x_continuous(name=NULL)

# p2a<-ggplot(tails,aes(x=tail_length,y=..density..,colour=sample))+
# geom_density(position='identity',size=.5)+theme_classic()+
# scale_y_continuous(name=NULL)+
# scale_x_continuous(name=NULL)

# p3<-ggplot(tails,aes(x=tail_length,y=..density..))+facet_wrap(~sample)+
# stat_bin(geom="step",binwidth=5,position='identity',size=.5)+theme_classic()+
# scale_y_continuous(name=NULL)+
# scale_x_continuous(name=NULL,limits=c(0,250))
# #theme(legend.position="none")

ggsave(plot=p1,file="/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/figures/version_9/Fig2Cpanel1.pdf",width=2,height=1.25) #PLOTS HERE
ggsave(plot=p2,file="/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/figures/version_9/Fig2Cpanel2.pdf",width=2,height=1.25) #PLOTS HERE



print("Calculating means ...")


colnames(tp40m_bs_mean)[1] <- "40 m"
colnames(tp1hr_bs_mean)[1] <- "1 h"
colnames(tp2hr_bs_mean)[1] <- "2 h"
colnames(tp4hr_bs_mean)[1] <- "4 h"
colnames(tp8hr_bs_mean)[1] <- "8 h"
colnames(tpSS_mean)[1]<-"SS"

m<-join_all(list(
tp40m_bs_mean,
tp1hr_bs_mean,
tp2hr_bs_mean,
tp4hr_bs_mean,
tp8hr_bs_mean,
tpSS_mean),by="accession",type="full")

m<-m[complete.cases(m),]
m_melt<-melt(m,id="accession")
# print(m[which(m$"40 m"<0),])
print(nrow(m))
m_melt[which(m_melt$value >= 250),]$value = 250
# m_melt[which(m_melt$value <= 0),]$value = 0

p3<-ggplot(m_melt,aes(x=value,y=..density..,colour=variable)) + 
	stat_bin(geom="step",binwidth=2,position='identity') + 
	scale_y_continuous(name=NULL,limits=c(0,NA),expand=c(0,0)) + 
	theme_tim() + 
	scale_x_continuous(name=NULL,limits=c(8,250), breaks = c(8,(1:5) * 50), labels = c(8,(1:5) * 50)) + 
	scale_colour_manual(values=cols) +
	theme(axis.line.x = element_blank()) + 
	geom_segment(aes(x = 8, xend = 250, y = 0, yend = 0))

ggsave(plot=p3,file="/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/figures/version_9/Fig2D.pdf",width=2,height=1.25) #PLOTS HERE