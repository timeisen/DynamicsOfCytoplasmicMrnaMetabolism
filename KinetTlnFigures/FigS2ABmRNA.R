#FigS2A and FigS2B
#Techincal and Biological Rep Comparisons
library(ggplot2)
library(reshape2)
library(plyr)
library(data.table)
source("ggplot_theme.R")
# theme_tim <- function(base_size = 6, base_family = "Helvetica")
# {
# theme_classic(base_size = base_size, base_family = base_family) +
# theme(
# 	panel.grid=element_blank(),
# 	axis.text.y = element_text (size = 6,color="black"),
# 	axis.text.x = element_text (size = 6,color="black"),
# 	axis.ticks.length = unit(0.1 , "cm"),
# 	#axis.line=element_line(size=.5,color="black"),
# 	#axis.text.x=element_blank(),
#   #axis.ticks=element_blank(),
#   #axis.title.x=element_blank(),
#   legend.position="none",
#   strip.background = element_blank(),
#   panel.spacing = unit(1, "lines"))
#   # axis.text.x=element_blank(),
#   # axis.text.y=element_blank(),
#   # axis.title.x=element_blank())
# }

stds <- c(as.character(read.table("/lab/solexa_bartel/teisen/RNAseq/Annotation_files/TAIL_seq/list_of_standards_TJE_notation.txt")$V1),
		  as.character(read.table("/lab/solexa_bartel/teisen/RNAseq/Annotation_files/TAIL_seq/list_of_standards_SWE_notation.txt")$V1),
		  as.character(read.table("/lab/solexa_bartel/teisen/RNAseq/Annotation_files/TAIL_seq/list_of_standards_AOS_notation.txt")$V1))



PAL_V1_1__<-fread("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/processed_files/miR-1_minus_aos_2h_means.txt",head=TRUE,sep="\t")
PAL_V1_155<-fread("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/processed_files/miR-155_minus_aos_2h_means.txt",head=TRUE,sep="\t")
PAL_V2_155<-fread("/lab/solexa_bartel/eichhorn/5EU_RPF/miR-1_miR-155_timecourse_tail_analysis/unfiltered_analysis/Tail_lengths/AATCCG-2_mediantails.txt",head=TRUE,sep="\t")
PAL_V2_1__<-fread("/lab/solexa_bartel/eichhorn/5EU_RPF/miR-1_miR-155_timecourse_tail_analysis/unfiltered_analysis/Tail_lengths/AATCCG-1_mediantails.txt",head=TRUE,sep="\t")
TAIL_155__<-fread("/lab/solexa_bartel/teisen/Tail-seq/miR-155_previous_analyses/miR-155_samples/2hr_minus_v8_pcr_17_03_03/CTAAGG-s_5_1_sequence_mediantails.txt",head=TRUE,sep="\t")
TAIL_1____<-fread("/lab/solexa_bartel/eichhorn/5EU_RPF/Old/miR-1_timecourse/170120/Tail_lengths_revised_4/GGCCAT-4_median_lengths.txt",head=TRUE,sep="\t")



colnames(PAL_V1_1__)[c(1,4)]<-c("accession","tags")
colnames(PAL_V1_155)[c(1,4)]<-c("accession","tags")
colnames(PAL_V2_155)[c(1,2)]<-c("accession","tags")
colnames(PAL_V2_1__)[c(1,2)]<-c("accession","tags")
colnames(TAIL_1____)[c(1,4)]<-c("accession","tags")
colnames(TAIL_155__)[c(1,4)]<-c("accession","tags")


PAL_V1_1__<-PAL_V1_1__[which(PAL_V1_1__$tags>=50),]
PAL_V1_155<-PAL_V1_155[which(PAL_V1_155$tags>=50),]
PAL_V2_155<-PAL_V2_155[which(PAL_V2_155$tags>=50),]
PAL_V2_1__<-PAL_V2_1__[which(PAL_V2_1__$tags>=50),]
TAIL_1____<-TAIL_1____[which(TAIL_1____$tags>=50),]
TAIL_155__<-TAIL_155__[which(TAIL_155__$tags>=50),]


PAL_V1_1__<-PAL_V1_1__[,c(1,3)]
PAL_V1_155<-PAL_V1_155[,c(1,3)]
PAL_V2_155<-PAL_V2_155[,c(1,4)]
PAL_V2_1__<-PAL_V2_1__[,c(1,4)]
TAIL_1____<-TAIL_1____[,c(1,2)]
TAIL_155__<-TAIL_155__[,c(1,2)]



colnames(PAL_V1_1__)[2]<-"PAL_V1_1__"
colnames(PAL_V1_155)[2]<-"PAL_V1_155"
colnames(PAL_V2_155)[2]<-"PAL_V2_155"
colnames(PAL_V2_1__)[2]<-"PAL_V2_1__"
colnames(TAIL_1____)[2]<-"TAIL_1____"
colnames(TAIL_155__)[2]<-"TAIL_155__"

m<-join_all(list(
	PAL_V1_1__, #2
	PAL_V1_155, #3
	PAL_V2_1__, #4
	PAL_V2_155, #5
	TAIL_1____, #6
	TAIL_155__),by="accession",type="full")  #7

m1 <- nrow(m[complete.cases(m[,c(2,3)]),])
m2 <- nrow(m[complete.cases(m[,c(6,7)]),])
m3 <- nrow(m[complete.cases(m[,c(4,5)]),])
m4 <- nrow(m[complete.cases(m[,c(4,6)]),])
m5 <- nrow(m[complete.cases(m[,c(2,4)]),])
m6 <- nrow(m[complete.cases(m[,c(2,6)]),])

print(m1)
print(m2)
print(m3)
print(m4)
print(m5)
print(m6)


print(nrow(m))
m <- m[which(!m$accession %in% stds),]
# m[m$PAL_V1_155>250 & !is.na(m$PAL_V1_155),]

print(cor(m[,-1],method="spearman",use="pairwise.complete.obs"))


p1<-ggplot(m,aes(x=PAL_V1_1__,y=PAL_V1_155))+geom_point(size=0.5,alpha=0.2)+theme_tim()+scale_x_continuous(limits=c(0,250),name=NULL,expand=c(0,0))+scale_y_continuous(limits=c(0,250),name=NULL,expand=c(0,0))+geom_abline(slope=1,intercept=0,linetype="dashed",color="grey")
p2<-ggplot(m,aes(x=TAIL_1____,y=TAIL_155__))+geom_point(size=0.5,alpha=0.2)+theme_tim()+scale_x_continuous(limits=c(0,250),name=NULL,expand=c(0,0))+scale_y_continuous(limits=c(0,250),name=NULL,expand=c(0,0))+geom_abline(slope=1,intercept=0,linetype="dashed",color="grey")
p3<-ggplot(m,aes(x=PAL_V2_1__,y=PAL_V2_155))+geom_point(size=0.5,alpha=0.2)+theme_tim()+scale_x_continuous(limits=c(0,250),name=NULL,expand=c(0,0))+scale_y_continuous(limits=c(0,250),name=NULL,expand=c(0,0))+geom_abline(slope=1,intercept=0,linetype="dashed",color="grey")
p4<-ggplot(m,aes(x=PAL_V2_1__,y=TAIL_1____))+geom_point(size=0.5,alpha=0.2)+theme_tim()+scale_x_continuous(limits=c(0,250),name=NULL,expand=c(0,0))+scale_y_continuous(limits=c(0,250),name=NULL,expand=c(0,0))+geom_abline(slope=1,intercept=0,linetype="dashed",color="grey")
p5<-ggplot(m,aes(x=PAL_V2_1__,y=PAL_V1_1__))+geom_point(size=0.5,alpha=0.2)+theme_tim()+scale_x_continuous(limits=c(0,250),name=NULL,expand=c(0,0))+scale_y_continuous(limits=c(0,250),name=NULL,expand=c(0,0))+geom_abline(slope=1,intercept=0,linetype="dashed",color="grey")
p6<-ggplot(m,aes(x=PAL_V1_1__,y=TAIL_1____))+geom_point(size=0.5,alpha=0.2)+theme_tim()+scale_x_continuous(limits=c(0,250),name=NULL,expand=c(0,0))+scale_y_continuous(limits=c(0,250),name=NULL,expand=c(0,0))+geom_abline(slope=1,intercept=0,linetype="dashed",color="grey")



ggsave(plot=p1,file="figures/version_7/FigS4A.pdf",width=2,height=2,useDingbats=FALSE)
ggsave(plot=p2,file="figures/version_7/FigS4B.pdf",width=2,height=2,useDingbats=FALSE)
ggsave(plot=p3,file="figures/version_7/FigS4C.pdf",width=2,height=2,useDingbats=FALSE)
ggsave(plot=p4,file="figures/version_7/FigS4D.pdf",width=2,height=2,useDingbats=FALSE)
ggsave(plot=p5,file="figures/version_7/FigS4E.pdf",width=2,height=2,useDingbats=FALSE)
ggsave(plot=p6,file="figures/version_7/FigS4F.pdf",width=2,height=2,useDingbats=FALSE)
