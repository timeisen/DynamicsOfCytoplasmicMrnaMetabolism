#fig1 with schwannhausser half lives, tail length comparison
library(tidyverse)
source("/lab/solexa_bartel/teisen/RNAseq/Scripts/general/ggplot_theme.R")

# hl<-read_tsv("/lab/solexa_bartel/teisen/RNAseq/Annotation_files/HL_schwan.txt")
# hl<-select(hl,`Refseq mRNA ID`,`mRNA half-life average [h]`)

hl <- read_tsv("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/processed_files/halflifeComparisons/miR-155_minu_PAL_halflives_logspace_global_offset_SS_RNAseqScaling.txt")

tails <- read_tsv("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/single_tag_files/background_subtracted_single_tag_files/SS_minus_norm_miR-155_v7_HYBRID20190731.txt",col_names = FALSE) #miR-155, means HYBRID
# tails <- read_tsv("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/single_tag_files/background_subtracted_single_tag_files/SS_minus_norm_miR-155_v7.txt",col_names = FALSE) #miR-155, means

IEGs<-read_tsv("/lab/solexa_bartel/eichhorn/5EU_RPF/miR-1_miR-155_analysis/Halflife_analysis/Halflife_measurements/miR-155_RNA_minusdox_alpha_hl_genenames_iegs.txt", col_names = FALSE)
RPGs<-read_tsv("processed_files/rpgs_list.txt",col_names = FALSE)
annot<-read_tsv("/lab/solexa_bartel/teisen/RNAseq/Annotation_files/TAIL_seq/symbol_to_accession.txt")

MeanTails <- tibble("accession" = tails$X252,
							 "mean" =  apply(tails[,-c(1:8,252)],1,function(x){weighted.mean(8:250, w = x)}))

colnames(hl)[1] <- "accession"

m<-inner_join(MeanTails,hl,by = "accession")
m <- inner_join(m, annot, by = "accession")

m$Type = "None"
colnames(IEGs)[c(1,2)] <- c("symbol","accession")
colnames(RPGs) <- "symbol"

m[which(m$accession %in% IEGs$accession),]$Type = "IEG"
m[which(m$symbol %in% RPGs$symbol),]$Type = "RPG"

m[m$mean < 50,]$mean = 50

m[m$halflife_t > 10^2,]$halflife_t = 10^2 
m[m$halflife_t < 10^-1,]$halflife_t = 10^-1
p1<-ggplot(m[which(m$Type == "None"),],aes(y = halflife_t,x = mean, color = Type))+
	geom_point(size = 0.5, alpha = 0.2,shape = 16)+
	geom_point(data = m[which(m$Type %in% c("IEG","RPG")),],aes(x = mean, y = halflife_t, color = Type),size = 0.5,shape = 16,alpha = 1) +
	scale_x_continuous(expand=c(0,0),trans = "log10",limits=c(50,250),breaks = c(50,100,150,200,250)) +
	scale_y_continuous(trans = "log10",limits = c(10^-1, 10^2),breaks = log_ticks(10^-1,10^2)[[2]],labels = log_ticks(10^-1,10^2)[[1]]) +
	scale_color_manual(values = c("red","black","lightskyblue")) +
	theme_tim() +
	theme(axis.line.y = element_blank()) + geom_segment(aes(x = 50, xend = 50, y = 0.1, yend = 1E2))

print(nrow(m))
print(cor(m$halflife_t,m$mean,method="spearman"))

ggsave(plot=p1,file="/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/figures/version_9/Fig1D.pdf",width=2,height=2,useDingbats=FALSE)