#FigS2F
#single_tails_minus_compile.R

library(ggplot2)
library(reshape2)
library(plyr)
binwidth = 1
bin_vector<-function(v,n){unname(tapply(v, (seq_along(v)-1) %/% n, sum))} #note that the 51st bin is the 251st tail length, it is excluded

source("/lab/solexa_bartel/teisen/RNAseq/Scripts/general/ggplot_theme.R")



library(ggplot2)
library(reshape2)
library(plyr)

stds<-read.table("/lab/solexa_bartel/teisen/RNAseq/Annotation_files/TAIL_seq/list_of_standards_SWE_notation.txt")

hl <- read.table("processed_files/halflifeComparisons/miR-155_minu_PAL_halflives_logspace_global_offset_SS.txt",head=TRUE)
colnames(hl)[5] <- "halflife"
# hl$halflife<-log(2)/hl$beta
# hl<-hl[,c(1,11)]
# print(quantile(hl$halflife))

raw<-read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/single_tag_files/40min_minus_tag_reformat_no_cutoff_miR-155.txt",sep="\t")
sub<-read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/single_tag_files/background_subtracted_single_tag_files/tp40m_bs_minus_miR-155_v7.txt",sep="\t")

#cleaning data frames
raw_accession<-raw[,1]
raw[,1]<-NULL
raw[,252]<-raw_accession

raw<-raw[which(rowSums(raw[,-252])>=50),]

colnames(raw)[252]<-"accession"
colnames(sub)[252]<-"accession"


long_hl_genes<-hl[hl$halflife>8,]$accession
short_hl_genes<-hl[hl$halflife<0.5,]$accession

raw_long_hl<-raw[which(raw$accession %in% long_hl_genes & raw$accession %in% sub$accession),]
raw_short_hl<-raw[which(raw$accession %in% short_hl_genes),]
sub_long_hl<-sub[which(sub$accession %in% long_hl_genes),]
sub_short_hl<-sub[which(sub$accession %in% short_hl_genes),]


print(nrow(raw_long_hl))
print(nrow(sub_long_hl))
print(nrow(raw_short_hl))
print(nrow(sub_short_hl))


tails<-data.frame(
	"raw_long_hl" = colSums(raw_long_hl[,-252]),
	"sub_long_hl" = colSums(sub_long_hl[,-252]),
	"raw_short_hl" = colSums(raw_short_hl[,-252]),
	"sub_short_hl" = colSums(sub_short_hl[,-252]))

tails<-data.frame(apply(tails,2,function(x){x/sum(x)}))
print(head(tails))
tails_long<-tails[,c(1,2)]
tails_short<-tails[,c(3,4)]
colnames(tails_long)[c(1,2)]<-c("raw","sub")
colnames(tails_short)[c(1,2)]<-c("raw","sub")

tails_long<-data.frame("raw" = bin_vector(tails_long[,1],binwidth),
						"sub" = bin_vector(tails_long[,2],binwidth))

tails_short<-data.frame("raw" = bin_vector(tails_short[,1],binwidth),
						"sub" = bin_vector(tails_short[,2],binwidth))

tails_long$hl<-"long"
tails_short$hl<-"short"
tails_long$tail_length<-0:(250/binwidth)
tails_short$tail_length<-0:(250/binwidth)


tails<-rbind(tails_long,tails_short)

tails_melt<-melt(tails,id=c("tail_length","hl"))


tails_melt$hl<-factor(tails_melt$hl,levels=c("short","long"))
tails_melt$size = "a"
tails_melt[which(tails_melt$variable == "sub"),]$size = "b"
tails_melt$variable <- factor(tails_melt$variable, levels = c("sub","raw"))
p1<-ggplot(tails_melt[which(tails_melt$tail_length*binwidth != 250),],aes(x=tail_length*binwidth,y=value,color=variable))+facet_wrap(~hl)+
geom_step(aes(size = size),size = rep(c(0.5,1),each = 500))+scale_y_continuous(name=NULL,expand=c(0,0),limits=c(-0.01,0.02))+theme_tim()+
scale_x_continuous(name=NULL,expand=c(0,0),limits=c(0,250))+geom_hline(yintercept=0,linetype="dashed",color="grey")

p2<-ggplot(tails,aes(x=raw,y=sub,))+facet_wrap(~hl)+
geom_point(size=0.5)+scale_y_continuous(trans="log10")+theme_tim()+
scale_x_continuous(trans="log10")+geom_abline(slope=1,intercept=0,linetype="dashed",color="grey")

ggsave(p1,file="figures/version_9/FigS4.pdf",width=2.5,height=2,useDingbats=FALSE)