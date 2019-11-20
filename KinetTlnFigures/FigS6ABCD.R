#halflife_figure
library(ggplot2)
library(plyr)
library(reshape2)
options("scipen"=5)
source("https://raw.githubusercontent.com/briatte/ggcorr/master/ggcorr.R")
source("/lab/solexa_bartel/teisen/RNAseq/Scripts/general/ggplot_theme.R")

#import 8 libraries
# miR_1xx_r0_minu<-read.table("processed_files/halflifeComparisons/miR-1_minu_R0_halflives_logspace_global_offset_noSS.txt",head=TRUE,sep="\t")
# miR_1xx_r0_plus<-read.table("processed_files/halflifeComparisons/miR-1_plus_R0_halflives_logspace_global_offset_noSS.txt",head=TRUE,sep="\t")
miR_1xx_pA_minu<-read.table("processed_files/halflifeComparisons/miR-1_minu_polyA_halflives_logspace_global_offset_SS.txt",head=TRUE,sep="\t")
miR_1xx_pA_plus<-read.table("processed_files/halflifeComparisons/miR-1_plus_polyA_halflives_logspace_global_offset_SS.txt",head=TRUE,sep="\t")
miR_1xx_PL_minu<-read.table("processed_files/halflifeComparisons/miR-1_minu_PAL_halflives_logspace_global_offset_SS_RNAseqScaling.txt",head=TRUE,sep="\t")
miR_1xx_PL_plus<-read.table("processed_files/halflifeComparisons/miR-1_plus_PAL_halflives_logspace_global_offset_SS_RNAseqScaling.txt",head=TRUE,sep="\t")
miR_155_pA_minu<-read.table("processed_files/halflifeComparisons/miR-155_minu_polyA_halflives_logspace_global_offset_SS.txt",head=TRUE,sep="\t")
miR_155_pA_plus<-read.table("processed_files/halflifeComparisons/miR-155_plus_polyA_halflives_logspace_global_offset_SS.txt",head=TRUE,sep="\t")
# miR_155_r0_minu<-read.table("processed_files/halflifeComparisons/miR-155_minu_R0_halflives_logspace_global_offset_noSS.txt",head=TRUE,sep="\t")
# miR_155_r0_plus<-read.table("processed_files/halflifeComparisons/miR-155_plus_R0_halflives_logspace_global_offset_noSS.txt",head=TRUE,sep="\t")
miR_155_PL_minu<-read.table("processed_files/halflifeComparisons/miR-155_minu_PAL_halflives_logspace_global_offset_SS_RNAseqScaling.txt",head=TRUE,sep="\t")
miR_155_PL_plus<-read.table("processed_files/halflifeComparisons/miR-155_plus_PAL_halflives_logspace_global_offset_SS_RNAseqScaling.txt",head=TRUE,sep="\t")

HL_schwan <- read.table("/lab/solexa_bartel/teisen/RNAseq/Annotation_files/HL_schwan_reformat.txt",head=TRUE,sep="\t")
HL_schwan <- HL_schwan[,c(1,4)]
colnames(HL_schwan) <- c("accession","Schwanhausser")


#get rid of genes that don't fit well
#??? Cutoff???
#COMPARE OFFSETS
print(head(miR_155_pA_minu))
# #trim
# miR_155_r0_minu<-miR_155_r0_minu[,c(1,5)]
# miR_155_r0_plus<-miR_155_r0_plus[,c(1,5)]
miR_155_pA_minu<-miR_155_pA_minu[,c(1,5)]
miR_155_pA_plus<-miR_155_pA_plus[,c(1,5)]
# miR_1xx_r0_minu<-miR_1xx_r0_minu[,c(1,5)]
# miR_1xx_r0_plus<-miR_1xx_r0_plus[,c(1,5)]
miR_1xx_pA_minu<-miR_1xx_pA_minu[,c(1,5)]
miR_1xx_pA_plus<-miR_1xx_pA_plus[,c(1,5)]
miR_155_PL_minu<-miR_155_PL_minu[,c(1,5)]
miR_155_PL_plus<-miR_155_PL_plus[,c(1,5)]
miR_1xx_PL_minu<-miR_1xx_PL_minu[,c(1,5)]
miR_1xx_PL_plus<-miR_1xx_PL_plus[,c(1,5)]

# colnames(miR_155_r0_minu)<-c("accession","miR_155_r0_minu")
# colnames(miR_155_r0_plus)<-c("accession","miR_155_r0_plus")
colnames(miR_155_pA_minu)<-c("accession","miR_155_pA_minu")
colnames(miR_155_pA_plus)<-c("accession","miR_155_pA_plus")
# colnames(miR_1xx_r0_minu)<-c("accession","miR_1xx_r0_minu")
# colnames(miR_1xx_r0_plus)<-c("accession","miR_1xx_r0_plus")
colnames(miR_1xx_pA_minu)<-c("accession","miR_1xx_pA_minu")
colnames(miR_1xx_pA_plus)<-c("accession","miR_1xx_pA_plus")
colnames(miR_155_PL_minu)<-c("accession","miR_155_PL_minu")
colnames(miR_155_PL_plus)<-c("accession","miR_155_PL_plus")
colnames(miR_1xx_PL_minu)<-c("accession","miR_1xx_PL_minu")
colnames(miR_1xx_PL_plus)<-c("accession","miR_1xx_PL_plus")


#Consider only no-site genes for the plus samples.
miR_1_nosite <- read.table("/archive/bartel/2015.03.04-11263/solexa_bartel/eichhorn/Ribosome_profiling_paper/Sites/Mouse/miR-1/mmu-miR-1_NoSite.txt")
miR_155_nosite <- read.table("/archive/bartel/2015.03.04-11263/solexa_bartel/eichhorn/Ribosome_profiling_paper/Sites/Mouse/miR-155/mmu-miR-155_NoSite.txt")

# miR_155_r0_plus <- miR_155_r0_plus[which(miR_155_r0_plus$accession %in% miR_155_nosite$V1),]
miR_155_pA_plus <- miR_155_pA_plus[which(miR_155_pA_plus$accession %in% miR_155_nosite$V1),]
miR_155_PL_plus <- miR_155_PL_plus[which(miR_155_PL_plus$accession %in% miR_155_nosite$V1),]
# miR_1xx_r0_plus <- miR_1xx_r0_plus[which(miR_1xx_r0_plus$accession %in% miR_1_nosite$V1),]
miR_1xx_pA_plus <- miR_1xx_pA_plus[which(miR_1xx_pA_plus$accession %in% miR_1_nosite$V1),]
miR_1xx_PL_plus <- miR_1xx_PL_plus[which(miR_1xx_PL_plus$accession %in% miR_1_nosite$V1),]

# print(nrow(miR_155_r0_minu))
# print(nrow(miR_155_r0_plus))
print(nrow(miR_155_pA_minu))
print(nrow(miR_155_pA_plus))
# print(nrow(miR_1xx_r0_minu))
# print(nrow(miR_1xx_r0_plus))
print(nrow(miR_1xx_pA_minu))
print(nrow(miR_1xx_pA_plus))
print(nrow(miR_155_PL_minu))
print(nrow(miR_155_PL_plus))
print(nrow(miR_1xx_PL_minu))
print(nrow(miR_1xx_PL_plus))

# print("miR_155_r0_minu")
# print(median(miR_155_r0_minu[,2]))
# print("miR_155_r0_plus")
# print(median(miR_155_r0_plus[,2]))
print("miR_155_pA_minu")
print(median(miR_155_pA_minu[,2]))
print("miR_155_pA_plus")
print(median(miR_155_pA_plus[,2]))
print("miR_1xx_r0_minu")
# print(median(miR_1xx_r0_minu[,2]))
# print("miR_1xx_r0_plus")
# print(median(miR_1xx_r0_plus[,2]))
# print("miR_1xx_pA_minu")
print(median(miR_1xx_pA_minu[,2]))
print("miR_1xx_pA_plus")
print(median(miR_1xx_pA_plus[,2]))
print("miR_155_PL_minu")
print(median(miR_155_PL_minu[,2]))
print("miR_155_PL_plus")
print(median(miR_155_PL_plus[,2]))
print("miR_1xx_PL_minu")
print(median(miR_1xx_PL_minu[,2]))
print("miR_1xx_PL_plus")
print(median(miR_1xx_PL_plus[,2]))

m<-join_all(list(
# miR_155_r0_minu,
# miR_155_r0_plus,
miR_155_pA_minu,
# miR_155_pA_plus,
# miR_1xx_r0_minu,
# miR_1xx_r0_plus,
miR_1xx_pA_minu,
# miR_1xx_pA_plus,
miR_155_PL_minu,
# miR_155_PL_plus,
miR_1xx_PL_minu,
# miR_1xx_PL_plus,
HL_schwan),type="full",by="accession")
print(nrow(m))
# m2<-m[complete.cases(m),]
m_melt<-melt(m,id="accession")

m[m$miR_155_pA_minu < 0.1 & !(is.na(m$miR_155_pA_minu < 0.1)),]$miR_155_pA_minu = 0.1
m[m$miR_1xx_pA_minu < 0.1 & !(is.na(m$miR_1xx_pA_minu < 0.1)),]$miR_1xx_pA_minu = 0.1
m[m$miR_155_PL_minu < 0.1 & !(is.na(m$miR_155_PL_minu < 0.1)),]$miR_155_PL_minu = 0.1
m[m$miR_1xx_PL_minu < 0.1 & !(is.na(m$miR_1xx_PL_minu < 0.1)),]$miR_1xx_PL_minu = 0.1

# m[m$miR_155_pA_minu > 1E2 & !(is.na(m$miR_155_pA_minu > 1E2)),]$miR_155_pA_minu = 1E2
# m[m$miR_1xx_pA_minu > 1E2 & !(is.na(m$miR_1xx_pA_minu > 1E2)),]$miR_1xx_pA_minu = 1E2
m[m$miR_155_PL_minu > 1E2 & !(is.na(m$miR_155_PL_minu > 1E2)),]$miR_155_PL_minu = 1E2
m[m$miR_1xx_PL_minu > 1E2 & !(is.na(m$miR_1xx_PL_minu > 1E2)),]$miR_1xx_PL_minu = 1E2


p3 <- ggcorr(data = NULL, cor_matrix = cor(log(m[, -1]), method="spearman",use="pairwise.complete.obs"),label=TRUE,label_round = 2,geom="text",hjust=0.75,layout.exp = 1,size = 2,label_size = 2)+theme_tim()
ggsave(plot=p3,file="/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/figures/version_9/FigS6C.pdf",width=4,height=4)

#corrplot(cor(m[,-1],method="spearman"),type="upper",method="number")
p1<-ggplot(m_melt[m_melt$variable %in% c("miR_155_r0_minu",
	"miR_155_r0_plus",
	"miR_155_pA_minu",
	"miR_155_pA_plus",
	"miR_155_PL_minu",
	"miR_155_PL_plus"),],aes(x=value,y=..density..,colour=variable))+stat_bin(geom="step",bins=50,position='identity')+theme_tim()+scale_x_continuous(trans="log10",name=NULL,expand=c(0,0),limits=c(0.001,100),labels=log_ticks(0.001,1E2)[[1]],breaks=log_ticks(0.001,1E2)[[2]])+scale_y_continuous(name=NULL,expand=c(0,0))+theme(legend.position="right")

p2<-ggplot(m_melt[m_melt$variable %in% c("miR_1xx_r0_minu",
	"miR_1xx_r0_plus",
	"miR_1xx_pA_minu",
	"miR_1xx_pA_plus",
	"miR_1xx_PL_minu",
	"miR_1xx_PL_plus"),],aes(x=value,y=..density..,colour=variable))+stat_bin(geom="step",bins=50,position='identity')+theme_tim()+scale_x_continuous(trans="log10",name=NULL,expand=c(0,0),limits=c(0.001,100),labels=log_ticks(0.001,1E2)[[1]],breaks=log_ticks(0.001,1E2)[[2]])+scale_y_continuous(name=NULL,expand=c(0,0))+theme(legend.position="right")

m3 <- merge(miR_155_PL_minu,HL_schwan,by="accession")




m3[m3$miR_155_PL_minu < 0.1,]$miR_155_PL_minu = 0.1
m3[m3$miR_155_PL_minu > 1E2,]$miR_155_PL_minu = 1E2
# m3[m3$Schwanhausser < 0.1,]$Schwanhausser = 0.1
# m3[m3$Schwanhausser > 1E2,]$Schwanhausser = 1E2



p4<-ggplot(m3,aes(x=miR_155_PL_minu,y=Schwanhausser))+geom_point(alpha=0.2,size=0.5,shape = 16)+theme_tim()+scale_x_continuous(name=NULL,trans="log10",limits=c(0.1,1E2),labels=log_ticks(0.1,1E2)[[1]],breaks=log_ticks(0.1,1E2)[[2]],expand=c(0,0))+scale_y_continuous(name=NULL,trans="log10",limits=c(0.1,1E2),labels=log_ticks(0.1,1E2)[[1]],breaks=log_ticks(0.1,1E2)[[2]],expand=c(0,0))+geom_abline(linetype="dashed",color="grey",slope=1,intercept=0)
print(cor(m3$miR_155_PL_minu,m3$Schwanhausser,method="spearman"))
print(nrow(m3))
print(head(m3))

ggsave(plot=p3,file="/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/figures/version_9/FigS6C.pdf",width=4,height=4)
ggsave(plot=p4,file="/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/figures/version_9/FigS6B2.pdf",width=2,height=2)
ggsave(plot=p1,width=4,height=2,file="/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/figures/version_9/FigS6B.pdf",useDingbats=FALSE)
ggsave(plot=p2,width=4,height=2,file="/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/figures/version_9/FigS6A.pdf",useDingbats=FALSE)