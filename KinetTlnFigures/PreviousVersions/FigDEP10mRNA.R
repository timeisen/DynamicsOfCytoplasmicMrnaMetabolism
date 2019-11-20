#DEP
library(ggplot2)
source("/lab/solexa_bartel/teisen/RNAseq/Scripts/general/ggplot_theme.R")


stds <- read.table("/lab/solexa_bartel/teisen/RNAseq/Annotation_files/TAIL_seq/list_of_standards_SWE_notation.txt")
expr<-read.table("/lab/solexa_bartel/teisen/RNAseq/kinetics_of_translation/actD_RNAseq/analysis/minus_counts_10rpm_cutoff_std3_normalized_0hr_cutoff_only.txt",head=TRUE)
expr <- expr[which(!expr$accession %in% stds$V1),]
expr[is.na(expr)] <- 0

m <- data.frame(
	time = c(0,1,3,7,15),
	expr = colSums(expr[,-1])/colSums(expr[,-1])[1])

p1 <- ggplot(m,aes(x=time,y=expr))+
geom_point(size=0.5)+
scale_x_continuous(expand=c(0,0))+
scale_y_continuous()+
theme_tim()

ggsave(plot=p1,file="figures/version_7/FigS7G.pdf",width=2,height=2,useDingbats=FALSE)