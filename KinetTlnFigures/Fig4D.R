
library(ggplot2)
options("scipen"=5)
source("/lab/solexa_bartel/teisen/RNAseq/Scripts/general/ggplot_theme.R")
options(warn=1)
hl <- read.table("processed_files/halflifeComparisons/miR-155_minu_PAL_halflives_logspace_global_offset_SS.txt",head=TRUE,sep="\t")
o8<-read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/rate_constant_measurements/miR-155_minus_samples_UNLINKV75_run6_final_global_param.txt")
annot<-read.table("/lab/solexa_bartel/teisen/RNAseq/Annotation_files/TAIL_seq/symbol_to_accession.txt",head=TRUE)
o8<-o8[complete.cases(o8),]


colnames(o8) <- c("accession","st_ul","a_ul","k_ul","b_ul","r_ul")
m<-merge(hl,o8,by="accession")
m<-merge(m,annot,by="accession")

print(nrow(m))
print(cor(m$k_ul,m$beta_t,method="spearman"))

m[which(m$halflife_t > 100),]$halflife_t = 100
print(m[which(m$k_ul > 1000 | m$k_ul <.01),])
p1<-ggplot(m,aes(x=halflife_t,y=k_ul,labels=symbol))+geom_point(size=0.5,alpha=.1)+
scale_x_continuous(name=NULL,trans="log10",limits=c(.001,100),breaks=log_ticks(.001,100)[[2]],labels=log_ticks(.001,100)[[1]],expand=c(.001,.001))+
scale_y_continuous(name=NULL,trans="log10",limits=c(.001,1000),breaks=log_ticks(.001,1000)[[2]],labels=log_ticks(.001,1000)[[1]],expand=c(.001,.001))+theme_tim()#expand_limits(x = .001, y = 0.0001)+


p1a<-ggplot(m,aes(x=halflife_t,y=k_ul,labels=symbol,color=-ratios))+geom_point()+
scale_x_continuous(name=NULL,trans="log10",limits=c(.001,100),breaks=log_ticks(.001,100)[[2]],labels=log_ticks(.001,100)[[1]],expand=c(.001,.001))+
scale_y_continuous(name=NULL,trans="log10",limits=c(.001,1000),breaks=log_ticks(.001,1000)[[2]],labels=log_ticks(.001,1000)[[1]],expand=c(.001,.001))#expand_limits(x = .001, y = 0.0001)+

m$ratios <- log2(m$halflife_t*m$k_ul)
m <- m[order(m$ratios),]
p2<-ggplot(m,aes(x=ratios))+geom_histogram(bins=100)+
scale_x_continuous()+
scale_y_continuous(expand=c(0,0))+
geom_vline(linetype="dashed",color="grey",xintercept=m[which(m$symbol %in% c("Pnrc1","1500012F01Rik","Hnrnph3")),]$ratios)+
theme_tim()

# ggsave(plot=p2,file="/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/figures/other/Fig4DHist.pdf",width=2,height=2)

break
ggsave(plot=p1,file="/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/figures/other/Fig4DwithSim.pdf",width=2,height=2)
