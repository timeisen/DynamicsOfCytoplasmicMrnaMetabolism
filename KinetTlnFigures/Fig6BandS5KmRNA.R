#Figure S5K, 6B
library(plyr)
library(reshape2)
library(tidyverse)
source("ggplot_theme.R")

dat <- read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/single_tag_files/background_subtracted_single_tag_files/SS_minus_norm_miR-155_v7_HYBRID20190731.txt")

hl<-read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/processed_files/halflifeComparisons/miR-155_minu_PAL_halflives_logspace_global_offset_SS_RNAseqScaling.txt",head=TRUE,sep="\t")
histones <- read.table("/lab/solexa_bartel/teisen/RNAseq/Annotation_files/hist_genes_mouse_accession_only.txt",head = FALSE)$V1
colnames(dat) <- c(250:0,"accession")
dat <- dat[which(!dat$accession %in% histones),]	
dat <- merge(dat,hl,by="accession")


long = 10
short = 20/60

ldf = dat[which(dat$halflife > long),2:252]
sdf = dat[which(dat$halflife < short),2:252]

print(nrow(ldf))
print(nrow(sdf))

# write.table(dat[which(dat$halflife < short),1],file = "/lab/solexa_bartel/teisen/Tail-seq/PalSeqDirectLig/Analysis/DirectTailCalls/short_hl_mRNA.txt",quote = FALSE, row.names = FALSE,sep = "\t",col.names = FALSE)
# write.table(dat[which(dat$halflife > long),1],file = "/lab/solexa_bartel/teisen/Tail-seq/PalSeqDirectLig/Analysis/DirectTailCalls/long_hl_mRNA.txt",quote = FALSE, row.names = FALSE,sep = "\t",col.names = FALSE)

# ldf <- ldf/rowSums(ldf) #uncomment for gene weighting.
# sdf <- sdf/rowSums(sdf) #uncomment for gene weighting.
ldfMeans <- apply(ldf,1,function(z){weighted.mean(z, x = 0:250, w = z)})
sdfMeans <- apply(sdf,1,function(z){weighted.mean(z, x = 0:250, w = z)})
MeansDistr <- tibble(
	accession = c(as.character(dat[which(dat$halflife > long),1]), 
		as.character(dat[which(dat$halflife < short),1])),
	mRNAmean = c(ldfMeans,sdfMeans),
	type = rep(c("Long half-life","Short half-life"),c(nrow(ldf),nrow(sdf)))
		)

ldfcs <- colSums(ldf)/sum(ldf)
sdfcs <- colSums(sdf)/sum(sdf)


print("Long half-life weighted mean:")
print(weighted.mean(x = 0:250, w = ldfcs))

print("Short half-life weighted mean:")
print(weighted.mean(x = 0:250, w = sdfcs))

dat155 <- data.frame(
	tail_length = rep(0:250,2),
	abundance = c(ldfcs, sdfcs),
	bin = rep(c("Long","Short"), each = 251))

#plotting

dat155$bin <- factor(dat155$bin, levels = c("Short","Long")) #plot short first
print(dat155[dat155$tail_length <= 50,])
print(dat155[dat155$tail_length == 250,])

p3 <- ggplot(dat155,aes(x = tail_length,y = abundance, color = bin)) + 
	  geom_step() +  
	  theme_tim() + 
	  #geom_smooth(data=dat,aes(x=tail_length,y=abundance),span=0.1,method="loess") + 
	  scale_y_continuous(limits=c(0,0.025),name=NULL,expand=c(0,0)) + 
	  scale_x_continuous(expand=c(0,0),limits=c(-0.5,250.5)) +
	  scale_color_manual(values = c("red","blue"))
p3b <- ggplot(dat155,aes(x = tail_length,y = abundance, color = bin)) + 
	  geom_step() +  
	  theme_tim() + 
	  #geom_smooth(data=dat,aes(x=tail_length,y=abundance),span=0.1,method="loess") + 
	  scale_y_continuous(limits=c(0,0.025),name=NULL,expand=c(0,0)) + 
	  scale_x_continuous(expand=c(0,0),limits=c(0,60)) +
	  scale_color_manual(values = c("red","blue"))

pV <- ggplot(MeansDistr, aes(x = type, y = mRNAmean, fill = type)) +
	geom_violin(draw_quantiles = c(0.5), width = 0.5, size = 0.5, color = "black") +
	# geom_jitter(height = 0, width = 0.2) + 
	scale_y_continuous(limits = c(0,200),expand = c(0,0)) + 
	scale_fill_manual(values = c("red","blue")) +
	theme_tim()

print("T test for violin plot:")
print(t.test(ldfMeans,sdfMeans,alternative = "two.sided")[["p.value"]])
print(paste("median sdf:", median(sdfMeans)))
print(paste("median ldf:", median(ldfMeans)))
ggsave(plot=pV,file="figures/version_9/FigS5Kviolin.pdf",width=2,height=2)
ggsave(plot=p3,file="figures/version_9/FigS8A.pdf",width = 40.0, height = 40.0, units = 'mm')
ggsave(plot=p3b,file="figures/version_9/FigS8AExpanded.pdf",width=2,height=2)



