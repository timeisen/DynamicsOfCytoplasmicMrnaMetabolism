#Mean tl vs hl comparisons. 2B and 7D
library(data.table)
library(plyr)
library(ggplot2)
library(reshape2)
source("/lab/solexa_bartel/teisen/RNAseq/Scripts/general/ggplot_theme.R")

#Read in the data
NMDcalls <- read.table("other_analyses/NMDsimulations/NMDcallsTailLengthAtDcp.txt",head=TRUE,sep="\t")
actD_1hr <- read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-1_actD_samples/analysis_20171030/model_input_files/global_normalized/actD_1hr_minus_background_subtracted_single_tag_v6.txt",sep="\t")
pulse_2hr <- read.table("single_tag_files/background_subtracted_single_tag_files/tp2hr_bs_minus_miR-155_v7.txt",sep="\t")
IEGs<-fread("/lab/solexa_bartel/eichhorn/5EU_RPF/miR-1_miR-155_analysis/Halflife_analysis/Halflife_measurements/miR-155_RNA_minusdox_alpha_hl_genenames_iegs.txt")
RPGs<-fread("processed_files/rpgs_list.txt",head=FALSE)
halflife<-fread("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/processed_files/halflifeComparisons/miR-155_minu_PAL_halflives_logspace_global_offset_SS_RNAseqScaling.txt",head=TRUE,sep="\t")
annot<-read.table("/lab/solexa_bartel/teisen/RNAseq/Annotation_files/TAIL_seq/symbol_to_accession.txt",head=TRUE)

#Clean the datasets
## Removed by TJE 2019 05 20
#  actD_1hr_t <-  actD_1hr[,-252]
# pulse_2hr_t <- pulse_2hr[,-252]
#  actD_1hr_t[ actD_1hr_t < 0] = 0
# pulse_2hr_t[pulse_2hr_t < 0] = 0
# print(dim(actD_1hr_t))

# print(dim(pulse_2hr_t))
#  actD_1hr[,-252] <-  actD_1hr_t
# pulse_2hr[,-252] <- pulse_2hr_t

actD_1hr$means <- apply(actD_1hr[,-252],1,function(z){weighted.mean(x=8:250,w=z[-c(1:8)])})
pulse_2hr$means <- apply(pulse_2hr[,-252],1,function(z){weighted.mean(x=8:250,w=z[-c(1:8)])})
actD_1hr <- actD_1hr[,c(252,253)]
pulse_2hr <- pulse_2hr[,c(252,253)]
IEGs<-IEGs[,2]
# RPGs<-RPGs[,2]
halflife <- halflife[,c(1,5)]

#Naming
colnames(actD_1hr) <- c("accession","MeanTl")
colnames(pulse_2hr) <- c("accession","MeanTl")
colnames(IEGs) <- "accession"
colnames(RPGs) <- "symbol"

#merge
actD_1hr <- merge(actD_1hr,halflife,by="accession")
pulse_2hr <- merge(pulse_2hr,halflife,by="accession")
actD_1hr <- merge(actD_1hr,annot,by="accession")
pulse_2hr <- merge(pulse_2hr,annot,by="accession")

#classification
actD_1hr$classification <- "None"
pulse_2hr$classification <- "None"
actD_1hr[which(actD_1hr$accession %in% IEGs$accession),]$classification = "IEG"
actD_1hr[which(actD_1hr$accession %in% NMDcalls$accession),]$classification = "NMD"
pulse_2hr[which(pulse_2hr$accession %in% IEGs$accession),]$classification = "IEG"
pulse_2hr[which(pulse_2hr$symbol %in% RPGs$symbol),]$classification = "RPG"
pulse_2hr[which(pulse_2hr$accession %in% NMDcalls$accession),]$classification = "NMD"


#Move points to axes
# pulse_2hr[pulse_2hr$MeanTl > 250,]$MeanTl = 250
pulse_2hr[pulse_2hr$MeanTl < 50,]$MeanTl = 50
# actD_1hr[actD_1hr$MeanTl > 250,]$MeanTl = 250

#plotting
pulse_2hr$classification <- factor(pulse_2hr$classification,levels = c("IEG","NMD","RPG","None"))
# pulse_2hr <- pulse_2hr[order(-pulse_2hr$classification),]


pulse_2hr[pulse_2hr$halflife_t > 10^2,]$halflife_t = 10^2 
pulse_2hr[pulse_2hr$halflife_t < 10^-1,]$halflife_t = 10^-1



p1 <- ggplot(pulse_2hr[which(pulse_2hr$classification %in% c("None","NMD")),],aes(x = MeanTl, y = halflife_t, color = classification,label = symbol))+
	geom_point(size = 0.5,shape = 16,alpha = 0.2)+
	geom_point(data = pulse_2hr[which(pulse_2hr$classification %in% c("IEG","RPG")),],aes(x = MeanTl, y = halflife_t, color = classification),size = 0.5,shape = 16,alpha = 1)+
	scale_x_continuous(expand=c(0,0),trans = "log10",limits=c(50,250),breaks = c(50,100,150,200,250))+
	scale_y_continuous(trans = "log10",limits = c(10^-1, 10^2),breaks = log_ticks(10^-1,10^2)[[2]],labels = log_ticks(10^-1,10^2)[[1]])+
	scale_color_manual(values = c("red","black","black","lightskyblue"))+
	theme_tim() +
	theme(axis.line.y = element_blank()) + geom_segment(aes(x = 50, xend = 50, y = 0.1, yend = 1E2))

actD_1hr[actD_1hr$halflife_t > 10^2,]$halflife_t = 10^2 
actD_1hr[actD_1hr$halflife_t < 10^-1,]$halflife_t = 10^-1 
actD_1hr[actD_1hr$MeanTl < 25,]$MeanTl = 25 

#Print correlations
print(cor(pulse_2hr$MeanTl,pulse_2hr$halflife_t,method="spearman"))
print(cor(actD_1hr$MeanTl,actD_1hr$halflife_t,method="spearman"))
print("Actd pearson cor")
print(cor(log(actD_1hr$MeanTl),log(actD_1hr$halflife_t),method="pearson"))
print(nrow(pulse_2hr))
print(nrow(actD_1hr))

p2 <- ggplot(actD_1hr,aes(x = MeanTl, y = halflife_t,label = symbol))+
	geom_point(size = 0.5, alpha = 0.2,shape = 16)+
	scale_x_continuous(expand=c(0,0),trans = "log10",limits=c(25,250),breaks = c(25,50,100,150,250))+
	scale_y_continuous(trans = "log10",limits = c(10^-1, 10^2),breaks = log_ticks(10^-1,10^2)[[2]],labels = log_ticks(10^-1,10^2)[[1]])+
	theme_tim() +
	theme(axis.line.y = element_blank()) + geom_segment(aes(x = 25, xend = 25, y = 0.1, yend = 1E2))


ggsave(plot = p1, file = "figures/version_9/Fig3BMeanTlHalflife.pdf",width = 2, height = 2, useDingbats = FALSE)
ggsave(plot = p2, file = "figures/version_9/Fig6BMeanTlHalflife.pdf",width = 2, height = 2, useDingbats = FALSE)

