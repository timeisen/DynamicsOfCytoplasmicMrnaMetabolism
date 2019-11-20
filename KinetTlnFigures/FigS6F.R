library(ggplot2)
library(data.table)
library(cowplot)
source("/lab/solexa_bartel/teisen/RNAseq/Scripts/general/ggplot_theme.R")

m <- fread("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/figures/other/SimVsDataLast7ntRmFromBothBoundedHYBRID20190326.txt")
colnames(m) <- c("accession","tp","data","sim")

# m[which(m$data>250),]$data = 250
# m[which(m$data<0),]$data = 0

p1 <- ggplot(m,aes(x=data,y=sim,color=as.factor(tp)))+
	geom_point(alpha=0.5,size=0.5,shape=16)+theme_tim()+
	scale_x_continuous(limits=c(0,250),expand=c(0,0))+
	scale_y_continuous(limits=c(0,250),expand=c(0,0),name=NULL)+
	geom_abline(color="grey",linetype="dashed")

legend <- get_legend(p1+theme(legend.position="left"))

m$diff <- m$data - m$sim
m$ratio <- log2(m$data/m$sim)

m <- m[order(-m$ratio),]
hl <- read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/processed_files/halflifeComparisons/miR-155_minu_PAL_halflives_logspace_global_offset_SS_RNAseqScaling.txt",head = TRUE)
m <- merge(m,hl)

p2 <- ggplot(m[which(m$tp == 5),],aes(x=data,y=sim,color=log10(halflife_t)))+
	geom_point(alpha=0.5,size=0.5,shape=16)+theme_tim()+
	scale_x_continuous(limits=c(0,250),expand=c(0,0))+
	scale_y_continuous(limits=c(0,250),expand=c(0,0),name=NULL)+
	geom_abline(color="grey",linetype="dashed") + 
	scale_color_gradient(low = "blue",high = "yellow")

legendP2 <- get_legend(p2+theme(legend.position="left"))

print("percentage of mRNAs with greater than 50 nt discrepancy at 40 min:")
print(nrow(m[which(m$tp == 5 & abs(m$diff) > 50),])/nrow(m[which(m$tp == 5),]))

ggsave(plot=p1,file="figures/version_8/SimVsDataMeanTL.pdf",width=2,height=2,useDingbats=FALSE)
ggsave(plot=legend,file="figures/version_8/LegendSimVsDataMeanTL.pdf",width=2,height=2,useDingbats=FALSE)
ggsave(plot=p2,file="figures/version_8/SimVsDataMeanTLFirstTpHalflife.pdf",width=2,height=2,useDingbats=FALSE)
ggsave(plot=legendP2,file="figures/version_8/SimVsDataMeanTLFirstTpHalflifelegend.pdf",width=2,height=2,useDingbats=FALSE)


print(cor(m$data,m$sim,method="spearman"))
print(cor(m$data,m$sim,method="pearson"))

print(length(unique(m$accession)))