# packrat::init("/lab/solexa_bartel/teisen/RNAseq/NeuroAnalysis/Analysis/notebooks/")
setwd("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/")
library(tidyverse)
source("/lab/solexa_bartel/teisen/RNAseq/Scripts/general/ggplot_theme.R")
library(ggpubr)

weightedSD = function(x){
    mn = weighted.mean(0:250,w = x)
    sdev = (sum(x * (0:250 - mn)^2)/sum(x))^0.5
    #returning only the sdev here
    return(sdev)}

hl <- read_tsv("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/processed_files/halflifeComparisons/miR-155_minu_PAL_halflives_logspace_global_offset_SS_RNAseqScaling.txt")
annot<-read.table("/lab/solexa_bartel/teisen/RNAseq/Annotation_files/TAIL_seq/symbol_to_accession.txt",head=TRUE)

#Note that there's a difference in the spearman corr when I use the hybrid datasets with the direct lig data. The direct lig data looks better, but
#both should be fine to use. 


#########UNCOMMENT THESE LINES TO RUN THE PAL-SEQ HYBRID DATA##########
###Read in and apply cutoffs to the data
miR1xx_uninduc_reformat <- read_tsv("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/single_tag_files/background_subtracted_single_tag_files/SS_minus_norm_miR-1_v7_HYBRID20190731.txt",col_names = FALSE)
miR155_uninduc_reformat <- read_tsv("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/single_tag_files/background_subtracted_single_tag_files/SS_minus_norm_miR-155_v7_HYBRID20190731.txt",col_names = FALSE)


#Format and process the data
colnames(miR1xx_uninduc_reformat) <- c(paste0("Tl_",0:250),"accession")
colnames(miR155_uninduc_reformat) <- c(paste0("Tl_",0:250),"accession")
miR1xx_uninduc_reformat <- miR1xx_uninduc_reformat[,c(252,1:251)]
miR155_uninduc_reformat <- miR155_uninduc_reformat[,c(252,1:251)]
#######################################################################

# ##########UNCOMMENT THESE LINES TO RUN THE DIRECT-LIG DATA#############
# miR1xx_uninduc_reformat <- read_delim("/lab/solexa_bartel/teisen/Tail-seq/PalSeqDirectLig/Data/Run3AllData/Tails20190603/miR1_uninduced_R0/tail_lengths_annot_end_only_CLEANED20190603_reformat.txt",delim="\t")
# miR155_uninduc_reformat <- read_delim("/lab/solexa_bartel/teisen/Tail-seq/PalSeqDirectLig/Data/Run3AllData/Tails20190603/miR155_uninduced_R0/tail_lengths_annot_end_only_CLEANED20190603_reformat.txt",delim="\t")
# miR1xx_uninduc_reformat <- miR1xx_uninduc_reformat[rowSums(miR1xx_uninduc_reformat[,-1])>=50,]
# miR155_uninduc_reformat <- miR155_uninduc_reformat[rowSums(miR155_uninduc_reformat[,-1])>=50,]
# ########################################################################

#Generate all cutoff data
miR155_uninduc_all_frac <- data.frame(frac = sapply(2:252,function(x){
    if(x==2){return(miR155_uninduc_reformat[,2]/rowSums(miR155_uninduc_reformat[,-1]))}
    rowSums(miR155_uninduc_reformat[,2:x])/rowSums(miR155_uninduc_reformat[,-1])}),
    accession = miR155_uninduc_reformat$accession)
colnames(miR155_uninduc_all_frac) <- c(paste0("frac",0:250),"accession")
miR155_uninduc_all_frac <- merge(miR155_uninduc_all_frac,hl[,c(1,5)],by="accession")

miR1xx_uninduc_all_frac <- data.frame(frac = sapply(2:252,function(x){
    if(x==2){return(miR1xx_uninduc_reformat[,2]/rowSums(miR1xx_uninduc_reformat[,-1]))}
    rowSums(miR1xx_uninduc_reformat[,2:x])/rowSums(miR1xx_uninduc_reformat[,-1])}),
    accession = miR1xx_uninduc_reformat$accession)
colnames(miR1xx_uninduc_all_frac) <- c(paste0("frac",0:250),"accession")
miR1xx_uninduc_all_frac <- merge(miR1xx_uninduc_all_frac,hl[,c(1,5)],by="accession")

miR1xx_uninduc_reformat_rev <- miR1xx_uninduc_reformat[,ncol(miR1xx_uninduc_reformat):1]
miR155_uninduc_reformat_rev <- miR155_uninduc_reformat[,ncol(miR155_uninduc_reformat):1]

#Generate all cutoff data, reversed
miR155_uninduc_all_frac_rev <- data.frame(frac = sapply(1:251,function(x){
    if(x==2){return(miR155_uninduc_reformat_rev[,1]/rowSums(miR155_uninduc_reformat_rev[,-252]))}
    rowSums(miR155_uninduc_reformat_rev[,1:x])/rowSums(miR155_uninduc_reformat_rev[,-252])}),
    accession = miR155_uninduc_reformat_rev$accession )
colnames(miR155_uninduc_all_frac_rev) <- c(paste0("frac",0:250),"accession")
miR155_uninduc_all_frac_rev <- merge(miR155_uninduc_all_frac_rev,hl[,c(1,5)],by="accession")

miR1xx_uninduc_all_frac_rev <- data.frame(frac = sapply(1:251,function(x){
    if(x==2){return(miR1xx_uninduc_reformat_rev[,1]/rowSums(miR1xx_uninduc_reformat_rev[,-252]))}
    rowSums(miR1xx_uninduc_reformat_rev[,1:x])/rowSums(miR1xx_uninduc_reformat_rev[,-252])}),
    accession = miR1xx_uninduc_reformat_rev$accession)
colnames(miR1xx_uninduc_all_frac_rev) <- c(paste0("frac",0:250),"accession")
miR1xx_uninduc_all_frac_rev <- merge(miR1xx_uninduc_all_frac_rev,hl[,c(1,5)],by="accession")



corFracVals <- tibble(
    tl = 0:250,
    tlrev = 250:0,
    miR1xx_cor = as.numeric(tail(cor(miR1xx_uninduc_all_frac[,-1],method="spearman",use="complete.obs"),1))[-252],
    miR155_cor = as.numeric(tail(cor(miR155_uninduc_all_frac[,-1],method="spearman",use="complete.obs"),1))[-252],
    miR1xx_cor_rev = as.numeric(tail(cor(miR1xx_uninduc_all_frac_rev[,-1],method="spearman",use="complete.obs"),1))[-252],
    miR155_cor_rev = as.numeric(tail(cor(miR155_uninduc_all_frac_rev[,-1],method="spearman",use="complete.obs"),1))[-252])



corFracValsGather <- gather(corFracVals,key = "DataType", value = "frac", -tl)
pCorFracVals <- ggplot(corFracValsGather[which(corFracValsGather$DataType %in% c("miR1xx_cor","miR155_cor")),], aes(x = tl, y = frac, color = DataType)) + 
    geom_step() + 
    scale_x_continuous(expand = c(0,0),limits = c(0,100),name = "Tail length cutoff") + 
    scale_y_continuous(expand = c(0,0),limits = c(-1,1),name = "Spearman R") + 
    theme_tim()

pCorFracValsRev <- ggplot(corFracValsGather[which(corFracValsGather$DataType %in% c("miR1xx_cor_rev","miR155_cor_rev")),], aes(x = tl, y = frac, color = DataType)) + 
    geom_step() + 
    scale_x_continuous(expand = c(0,0),limits = c(0,150),name = "Tail length cutoff",breaks = (0:6)*25,labels = 250 - (0:6)*25) + 
    scale_y_continuous(expand = c(0,0),limits = c(-1,1),name = "Spearman R") + 
    theme_tim()


miR155_uninduc_short_tail_frac <- miR155_uninduc_reformat
miR155_uninduc_short_tail_frac$frac <- rowSums(miR155_uninduc_short_tail_frac[,2:20])/rowSums(miR155_uninduc_short_tail_frac[,-1])
miR155_uninduc_short_tail_frac <- subset(miR155_uninduc_short_tail_frac,select = c("accession","frac"))

miR155_uninduc_short_tail_frac_symbol <- inner_join(miR155_uninduc_short_tail_frac,annot,by = "accession")
miR155_uninduc_short_tail_frac_symbol <- miR155_uninduc_short_tail_frac_symbol[order(miR155_uninduc_short_tail_frac_symbol$frac,decreasing = TRUE),]
print(miR155_uninduc_short_tail_frac_symbol,n = 30)

miR1xx_uninduc_short_tail_frac <- miR1xx_uninduc_reformat
miR1xx_uninduc_short_tail_frac$frac <- rowSums(miR1xx_uninduc_short_tail_frac[,2:20])/rowSums(miR1xx_uninduc_short_tail_frac[,-1])
miR1xx_uninduc_short_tail_frac <- subset(miR1xx_uninduc_short_tail_frac,select = c("accession","frac"))

miR155_uninduc_sd <- miR155_uninduc_reformat
miR155_uninduc_sd$sd <- apply(miR155_uninduc_sd[,-1],1,weightedSD)
miR155_uninduc_sd <- subset(miR155_uninduc_sd,select = c("accession","sd"))


miR155_uninduc_short_tail_frac_rev <- miR155_uninduc_reformat_rev
miR155_uninduc_short_tail_frac_rev$frac <- rowSums(miR155_uninduc_short_tail_frac_rev[,1:75])/rowSums(miR155_uninduc_short_tail_frac_rev[,-252])
miR155_uninduc_short_tail_frac_rev <- subset(miR155_uninduc_short_tail_frac_rev,select = c("accession","frac"))

miR155_uninduc_frac_hl <- inner_join(miR155_uninduc_short_tail_frac,hl,by="accession")

miR1xx_uninduc_frac_hl <- inner_join(miR1xx_uninduc_short_tail_frac,hl,by="accession")
# print(cor(miR1xx_uninduc_frac_hl$frac,miR1xx_uninduc_frac_hl$halflife_t,method="spearman",use = "complete.obs"))
# print(dim(miR1xx_uninduc_frac_hl))

miR155_uninduc_short_tail_frac_rev_hl <- inner_join(miR155_uninduc_short_tail_frac_rev,hl,by="accession")


miR155_uninduc_sd_hl <- inner_join(miR155_uninduc_sd,hl,by="accession")


miR155_uninduc_frac_hl[miR155_uninduc_frac_hl$halflife_t < 0.1,]$halflife_t = 0.1
miR155_uninduc_frac_hl[miR155_uninduc_frac_hl$halflife_t > 100,]$halflife_t = 100

miR155_uninduc_frac_hl[miR155_uninduc_frac_hl$frac < 0.001,]$frac = 0.001

pmiR155_frac <- ggplot(miR155_uninduc_frac_hl,aes(x = halflife_t, y = frac)) + 
    geom_point(alpha = 0.2,size = 0.5, shape = 16) + 
    scale_x_continuous(trans = "log10", limits = c(0.1,100),breaks=log_ticks(0.1,100)[[2]],labels=log_ticks(0.1,100)[[1]],name = "Half-life") + 
    scale_y_continuous(trans = "log10", limits = c(0.001,1),breaks=log_ticks(0.001,1)[[2]],labels=log_ticks(0.001,1)[[1]],name = "Fraction < 20") +
    # ggtitle("Short tailed miR-155 uninduced mRNAs and half-life") +
    # stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top",size = 3, aes(label = paste(..r.label..)))+
    theme_tim() +
    theme(axis.line.x = element_blank()) + geom_segment(aes(x = 0.1, xend = 100, y = 0, yend = 0)) +
    theme(axis.line.y = element_blank()) + 
    geom_segment(aes(y = 0.001, yend = 1, x = 0, xend = 0))

print("fraction comparison n:")
pmiR155_frac_histogram <- ggplot(miR155_uninduc_short_tail_frac,aes(x = frac, y = ..count../sum(..count..))) + 
    stat_bin(geom="step",bins=100,position='identity') +
    scale_x_continuous(expand=c(0,0)) + 
    scale_y_continuous(expand=c(0,0)) +
    # ggtitle("Short tailed miR-155 uninduced mRNAs and half-life") +
    # stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top",size = 3, aes(label = paste(..r.label..)))+
    theme_tim()

pmiR155_frac_CDF <- ggplot(miR155_uninduc_short_tail_frac,aes(x = frac)) + 
    stat_ecdf() +
    scale_x_continuous(expand=c(0,0),breaks = c(0:4)/4,limits = c(0,1)) + 
    scale_y_continuous(expand=c(0,0),breaks = c(0:4)/4) +
    # ggtitle("Short tailed miR-155 uninduced mRNAs and half-life") +
    # stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top",size = 3, aes(label = paste(..r.label..)))+
    theme_tim()

pmiR1xx_frac <- ggplot(miR1xx_uninduc_frac_hl,aes(x = halflife_t, y = frac)) + 
    geom_point(alpha = 0.2,size = 0.5, shape = 16) + 
    scale_x_continuous(trans = "log10", limits = c(0.01,100),breaks=log_ticks(0.01,100)[[2]],labels=log_ticks(0.01,100)[[1]],expand=c(0,0),name = "Half-life") + 
    scale_y_continuous(trans = "log10", limits = c(0.001,1),breaks=log_ticks(0.001,1)[[2]],labels=log_ticks(0.001,1)[[1]],expand=c(0,0),name = "Fraction < 20") +
    # ggtitle("Short tailed miR-1 uninduced mRNAs and half-life") +
    # stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top",size = 3, aes(label = paste(..r.label..)))+
    theme_tim() +
    theme(axis.line.x = element_blank()) + geom_segment(aes(x = 0.1, xend = 100, y = 0, yend = 0))

#175mer
miR155_uninduc_short_tail_frac_rev_hl[miR155_uninduc_short_tail_frac_rev_hl$halflife_t < 0.1,]$halflife_t = 0.1
miR155_uninduc_short_tail_frac_rev_hl[miR155_uninduc_short_tail_frac_rev_hl$halflife_t > 100,]$halflife_t = 100
pmiR155_frac_rev <- ggplot(miR155_uninduc_short_tail_frac_rev_hl,aes(x = halflife_t, y = frac)) + 
    geom_point(alpha = 0.2,size = 0.5, shape = 16) + 
    scale_x_continuous(trans = "log10", limits = c(0.1,100),breaks=log_ticks(0.1,100)[[2]],labels=log_ticks(0.1,100)[[1]],name = "Half-life") + 
    scale_y_continuous(trans = "log10", limits = c(0.001,1),breaks=log_ticks(0.001,1)[[2]],labels=log_ticks(0.001,1)[[1]],name = "Fraction < 20") +
    # ggtitle("Short tailed miR-155 uninduced mRNAs and half-life") +
    # stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top",size = 3, aes(label = paste(..r.label..)))+
    theme_tim() +
    theme(axis.line.x = element_blank()) + geom_segment(aes(x = 0.1, xend = 100, y = 0, yend = 0)) +
    theme(axis.line.y = element_blank()) + 
    geom_segment(aes(y = 0.001, yend = 1, x = 0, xend = 0))

miR155_uninduc_sd_hl[miR155_uninduc_sd_hl$halflife_t < 0.1,]$halflife_t = 0.1
miR155_uninduc_sd_hl[miR155_uninduc_sd_hl$halflife_t > 100,]$halflife_t = 100

# miR155_uninduc_sd_hl[miR155_uninduc_sd_hl$sd < 25,]$sd = 25
miR155_uninduc_sd_hl[miR155_uninduc_sd_hl$sd > 75,]$sd = 75

pmiR155_sd <- ggplot(miR155_uninduc_sd_hl,aes(x = halflife_t, y = sd)) + 
    geom_point(alpha = 0.2,size = 0.5, shape = 16) + 
    scale_x_continuous(trans = "log10", limits = c(0.1,100),breaks=log_ticks(0.1,100)[[2]],labels=log_ticks(0.1,100)[[1]],name = "Half-life") + 
    scale_y_continuous(limits = c(25,75),breaks = c(25,40,60,75),name = "SD") +
    # ggtitle("Short tailed miR-155 uninduced mRNAs and half-life") +
    # stat_cor(method = "spearman", label.x.npc = "left", label.y.npc = "top",size = 3, aes(label = paste(..r.label..)))+
    theme_tim() +
    theme(axis.line.y = element_blank()) + 
    # geom_segment(aes(x = 0.1, xend = 100, y = 25, yend = 25)) +
    geom_segment(aes(y = 25, yend = 75, x = 0, xend = 0))

print("short tail frac vs halflife comparison:")
print(cor(miR155_uninduc_frac_hl$frac,miR155_uninduc_frac_hl$halflife_t,method="spearman",use = "complete.obs"))
print(dim(miR155_uninduc_frac_hl))
print("sd vs halflife comparison:")
print(cor(miR155_uninduc_sd_hl$sd,miR155_uninduc_sd_hl$halflife_t,method="spearman",use = "complete.obs"))
print(dim(miR155_uninduc_sd_hl))
print("long tail frac vs halflife comparison:")
print(cor(miR155_uninduc_short_tail_frac_rev_hl$frac,miR155_uninduc_short_tail_frac_rev_hl$halflife_t,method="spearman",use = "complete.obs"))
print(dim(miR155_uninduc_short_tail_frac_rev_hl))

ggsave(plot = pmiR155_sd,   file = "figures/version_9/pCorFracValsFigS11F.pdf", units = 'mm', width = 40, height = 40, useDingbats = FALSE)
ggsave(plot = pCorFracVals, file = "figures/version_9/pCorFracValsFigS11A.pdf", units = 'mm', width = 40, height = 40, useDingbats = FALSE)
ggsave(plot = pmiR155_frac, file = "figures/version_9/pCorFracValsFigS11B.pdf", units = 'mm', width = 40, height = 40, useDingbats = FALSE)
# ggsave(plot = pmiR1xx_frac, file = "figures/version_9/pCorFracValsFigS11C.pdf", units ='mm', width = 40, height = 40, useDingbats = FALSE)
ggsave(plot = pmiR155_frac_rev, file = "figures/version_9/pCorFracValsFigS11D.pdf", units = 'mm', width = 40, height = 40, useDingbats = FALSE)
ggsave(plot = pCorFracValsRev,  file = "figures/version_9/pCorFracValsFigS11E.pdf", units = 'mm', width = 40, height = 40, useDingbats = FALSE)
ggsave(plot = pmiR155_frac_histogram,  file = "figures/version_9/pCorFracValsHistogram.pdf", units = 'mm', width = 40, height = 40, useDingbats = FALSE)
ggsave(plot = pmiR155_frac_CDF,  file = "figures/version_9/pCorFracValsCDF.pdf", units = 'mm', width = 40, height = 40, useDingbats = FALSE)