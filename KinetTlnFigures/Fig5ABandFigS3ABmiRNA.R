#Fig5AB and FigS3AB miRNA

library(ggplot2)
source("/lab/solexa_bartel/teisen/RNAseq/Scripts/general/ggplot_theme.R")
# minu <- read.table("/lab/solexa_bartel/teisen/RNAseq/kinetics_of_translation/final_analyses/RNA/miR-155/Total-RNA_17_06_17/CAGATC-s_2_1_expression_values_10rpm_cutoff.txt",head=TRUE)
# plus <- read.table("/lab/solexa_bartel/teisen/RNAseq/kinetics_of_translation/final_analyses/RNA/miR-155/Total+RNA_17_06_17/AGCGCT-s_8_1_expression_values.txt",head=TRUE)

# miR-155_rates_fc_adjusted.txt:
#         python /lab/bartel4_ata/Stephen/correct_utrlen_effect.py /lab/solexa_bartel/eichhorn/5EU_RPF/miR-1_miR-155_analysis/Corrected_foldchanges/mm10_tpUTR_len.txt /archive/bartel/2015.03.04-11263/solexa_bartel/eichhorn/Ribosome_profiling_paper/Sites/Mouse/miR-155/mmu-miR-155_NoSite.txt SSfcValues.txt SSfcValuesAdjusted.txt

# miR-155_nosite_rates.txt miR-155_tputr_rates.txt miR-155_topsite_rates.txt: miR-155_rates_fc_adjusted.txt
#         python /lab/bartel4_ata/Stephen/divide_RPKM_by_sites.py 3 /archive/bartel/2015.03.04-11263/solexa_bartel/eichhorn/Ribosome_profiling_paper/Sites/Mouse/miR-155/mmu-miR-155_NoSite.txt /archive/bartel/2015.03.04-11263/solexa_bartel/eichhorn/Ribosome_profiling_paper/Sites/Mouse/miR-155/mmu-miR-155_tpUTR_oneormore.txt /lab/solexa_bartel/eichhorn/5EU_paper/Eichhorn_data/RESTORE/miR-155_3T3_repressed_genes_25pct_down.txt miR-155_nosite_SSfc.txt miR-155_tputr_SSfc.txt miR-155_topsite_SSfc.txt SSfcValuesAdjusted.txt

## miR-1
# minu <- read.table("/lab/solexa_bartel/teisen/RNAseq/kinetics_of_translation/final_analyses/RNA/miR-1/Total-RNA_17_06_17/CAGATC-s_2_1_expression_values_10rpm_cutoff.txt",head=TRUE)
# plus <- read.table("/lab/solexa_bartel/teisen/RNAseq/kinetics_of_translation/final_analyses/RNA/miR-1/Total+RNA_17_06_17/CAGATC-s_3_1_expression_values.txt",head=TRUE)

# miR-1_rates_fc_adjusted.txt: miR-1_rates_fc.txt
# 	python /lab/bartel4_ata/Stephen/correct_utrlen_effect.py /lab/solexa_bartel/eichhorn/5EU_RPF/miR-1_miR-155_analysis/Corrected_foldchanges/mm10_tpUTR_len.txt /archive/bartel/2015.03.04-11263/solexa_bartel/eichhorn/Ribosome_profiling_paper/Sites/Mouse/miR-1/mmu-miR-1_NoSite.txt miR1_SSfcValues.txt miR1_SSfcValuesAdjusted.txt

# miR-1_nosite_rates.txt miR-1_tputr_rates.txt miR-1_topsite_rates.txt: miR-1_rates_fc_adjusted.txt
# 	python /lab/bartel4_ata/Stephen/divide_RPKM_by_sites.py 3 /archive/bartel/2015.03.04-11263/solexa_bartel/eichhorn/Ribosome_profiling_paper/Sites/Mouse/miR-1/mmu-miR-1_NoSite.txt /archive/bartel/2015.03.04-11263/solexa_bartel/eichhorn/Ribosome_profiling_paper/Sites/Mouse/miR-1/mmu-miR-1_tpUTR_oneormore.txt /lab/solexa_bartel/eichhorn/5EU_paper/Eichhorn_data/miR-1_3T3_repressed_genes_25pct_down.txt miR-1_nosite.txt miR-1_tputr.txt miR-1_topsite.txt miR1_SSfcValuesAdjusted.txt

allFCmiR155 <- read.table("SSfcValuesAdjusted.txt",head=TRUE)
noSitemiR155 <- read.table("miR-155_nosite_SSfc.txt",head=TRUE)
# topSitemiR155 <- read.table("miR-155_topsite_SSfc.txt",head=TRUE)
sitemiR155 <- read.table("miR-155_tputr_SSfc.txt",head=TRUE)

s155pp <- read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/other_analyses/directLig/PalAnnot3pEndsmiR/TargetScan7.1__miR-155-5p.predicted_targets.txt",head=TRUE,sep="\t")
print(nrow(s155pp))
colnames(s155pp)[1] <- "symbol"
annot<-read.table("/lab/solexa_bartel/teisen/RNAseq/Annotation_files/TAIL_seq/symbol_to_accession.txt",head=TRUE)
sitemiR155 <- merge(sitemiR155,annot)
print(nrow(sitemiR155))

print(sitemiR155[!sitemiR155$symbol %in% s155pp$symbol,]$symbol)

sitemiR155 <- merge(sitemiR155,s155pp)
print(nrow(sitemiR155))
colnames(sitemiR155)[5] = "TotalCppScore"

allFCmiR155$site <- "none"
#Order matters for the next 3 lines. 
allFCmiR155[which(allFCmiR155$accession %in% noSitemiR155$accession),]$site = "NoSite"
allFCmiR155[which(allFCmiR155$accession %in% sitemiR155$accession),]$site = "Site"
# allFCmiR155[which(allFCmiR155$accession %in% topSitemiR155$accession),]$site = "TopSite"


hl <- read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/processed_files/halflifeComparisons/miR-155_minu_PAL_halflives_logspace_global_offset_SS.txt",head=TRUE)
hl <- hl[,c(1,5)]

allFCmiR155 <- merge(allFCmiR155,hl,by="accession")

#Subtract the median of the NoSite
allFCmiR155[allFCmiR155$site == "NoSite",]$fcSS = allFCmiR155[allFCmiR155$site == "NoSite",]$fcSS - median(allFCmiR155[allFCmiR155$site == "NoSite",]$fcSS)
allFCmiR155[allFCmiR155$site == "Site",]$fcSS = allFCmiR155[allFCmiR155$site == "Site",]$fcSS - median(allFCmiR155[allFCmiR155$site == "NoSite",]$fcSS)
# allFCmiR155[allFCmiR155$site == "TopSite",]$fcSS = allFCmiR155[allFCmiR155$site == "TopSite",]$fcSS - median(allFCmiR155[allFCmiR155$site == "NoSite",]$fcSS)

allFCmiR155 <- allFCmiR155[order(allFCmiR155$site),]
# allFCmiR155 <- allFCmiR155[which(allFCmiR155$site %in% c("Site","TopSite")),]

NoSitemiR155 <- allFCmiR155[allFCmiR155$site == "NoSite",]
SitemiR155 <- allFCmiR155[allFCmiR155$site == "Site",]
# TopSitemiR155 <- allFCmiR155[allFCmiR155$site == "TopSite",]

# print(cor(allFCmiR155[which(!allFCmiR155$site=="none"),]$fcSS,allFCmiR155[which(!allFCmiR155$site=="none"),]$halflife_t,method="spearman"))
# print(nrow(allFCmiR155[which(!allFCmiR155$site=="none"),]))
print("NO SITE miR-155")
print(nrow(NoSitemiR155))
print(cor.test(SitemiR155$halflife_t,SitemiR155$fcSS,method="spearman"))
print(nrow(SitemiR155))
print(t.test(NoSitemiR155$halflife_t,SitemiR155$halflife_t,alternative="greater"))
# print(cor.test(TopSitemiR155$halflife_t,TopSitemiR155$fcSS,method="spearman"))
# print(nrow(TopSitemiR155))
# print(t.test(NoSitemiR155$halflife_t,TopSitemiR155$halflife_t,alternative="greater"))


allFCmiR155[which(allFCmiR155$accession=="NM_008597"),]$fcSS = -2
p1 <- ggplot(allFCmiR155[which(!allFCmiR155$site=="none"),],aes(x=halflife_t,y=fcSS,color=site))+
	geom_point(size=0.5)+
	scale_x_continuous(trans="log2",expand=c(0,0))+
	theme_tim()+
	scale_y_continuous(limits=c(-2.0,2.0),name=NULL,expand=c(0,0))+
	geom_hline(yintercept=median(allFCmiR155[which(allFCmiR155$site %in% c("Site")),]$fcSS),color="green")+
	# geom_hline(yintercept=median(allFCmiR155[which(allFCmiR155$site %in% c("TopSite")),]$fcSS),color="blue")+
	geom_hline(yintercept=median(allFCmiR155[which(allFCmiR155$site %in% c("NoSite")),]$fcSS),color="red")+
	geom_vline(xintercept=median(allFCmiR155[which(allFCmiR155$site %in% c("Site")),]$halflife_t),color="green")+
	# geom_vline(xintercept=median(allFCmiR155[which(allFCmiR155$site %in% c("TopSite")),]$halflife_t),color="blue")+
	geom_vline(xintercept=median(allFCmiR155[which(allFCmiR155$site %in% c("NoSite")),]$halflife_t),color="red")

# ggsave(plot=p1,file="/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/figures/version_6/FigS7X.pdf",width=2,height=2,useDingbats=FALSE)

##miR-1

allFCmiR1 <- read.table("miR1_SSfcValuesAdjusted.txt",head=TRUE)
noSitemiR1 <- read.table("miR-1_nosite_SSfc.txt",head=TRUE)
# topSitemiR1 <- read.table("miR-1_topsite_SSfc.txt",head=TRUE)
sitemiR1 <- read.table("miR-1_tputr_SSfc.txt",head=TRUE)

s1pp <- read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/other_analyses/directLig/PalAnnot3pEndsmiR/TargetScan7.1__miR-1-3p_206-3p.predicted_targets.txt",head=TRUE,sep="\t")
print(nrow(s1pp))
colnames(s1pp)[1] <- "symbol"
sitemiR1 <- merge(sitemiR1,annot)
# print(nrow(sitemiR1))
sitemiR1 <- merge(sitemiR1,s1pp)
# print(nrow(sitemiR1))
colnames(sitemiR1)[5] = "TotalCppScore"
sitemiR155 <- merge(sitemiR155,hl)
sitemiR1   <- merge(sitemiR1,hl)
print("Numbers of sites, context ++")
sitemiR1$TotalCppScore <- as.numeric(as.character(sitemiR1$TotalCppScore))
sitemiR155$TotalCppScore <- as.numeric(as.character(sitemiR155$TotalCppScore))



allFCmiR1$site <- "none"
#Order matters for the next 3 lines. 
allFCmiR1[which(allFCmiR1$accession %in% noSitemiR1$accession),]$site = "NoSite"
allFCmiR1[which(allFCmiR1$accession %in% sitemiR1$accession),]$site = "Site"
# allFCmiR1[which(allFCmiR1$accession %in% topSitemiR1$accession),]$site = "TopSite"


allFCmiR1 <- merge(allFCmiR1,hl,by="accession")

#Subtract the median of the NoSite
allFCmiR1[allFCmiR1$site == "NoSite",]$ssFC = allFCmiR1[allFCmiR1$site == "NoSite",]$ssFC - median(allFCmiR1[allFCmiR1$site == "NoSite",]$ssFC)
allFCmiR1[allFCmiR1$site == "Site",]$ssFC = allFCmiR1[allFCmiR1$site == "Site",]$ssFC - median(allFCmiR1[allFCmiR1$site == "NoSite",]$ssFC)
# allFCmiR1[allFCmiR1$site == "TopSite",]$ssFC = allFCmiR1[allFCmiR1$site == "TopSite",]$ssFC - median(allFCmiR1[allFCmiR1$site == "NoSite",]$ssFC)

allFCmiR1 <- allFCmiR1[order(allFCmiR1$site),]
# allFCmiR1 <- allFCmiR1[which(allFCmiR1$site %in% c("Site","TopSite")),]

NoSitemiR1 <- allFCmiR1[allFCmiR1$site == "NoSite",]
SitemiR1 <- allFCmiR1[allFCmiR1$site == "Site",]
# TopSitemiR1 <- allFCmiR1[allFCmiR1$site == "TopSite",]

# print(cor(allFCmiR1[which(!allFCmiR1$site=="none"),]$ssFC,allFCmiR1[which(!allFCmiR1$site=="none"),]$halflife_t,method="spearman"))
# print(nrow(allFCmiR1[which(!allFCmiR1$site=="none"),]))
print("NO SITE miR-1")
print(nrow(NoSitemiR1))

print(cor.test(SitemiR1$halflife_t,SitemiR1$ssFC,method="spearman"))
print(nrow(SitemiR1))
# print(cor.test(TopSitemiR1$halflife_t,TopSitemiR1$ssFC,method="spearman"))
# print(nrow(TopSitemiR1))


p2 <- ggplot(allFCmiR1[which(!allFCmiR1$site %in% c("none","NoSite")),],aes(x=halflife_t,y=ssFC,color=site))+geom_point()+
	scale_x_continuous(trans="log2",expand=c(0,0),limits=2^c(-9.08514,6.97736),breaks=2^c(-6,-2,2,6),labels=c("0.015625","0.25","4","64"))+theme_tim()+
	scale_y_continuous(limits=c(-2.0,2.0),name=NULL,expand=c(0,0))+
	geom_hline(yintercept=median(allFCmiR1[which(allFCmiR1$site %in% c("Site")),]$ssFC),color="green")+
	# geom_hline(yintercept=median(allFCmiR1[which(allFCmiR1$site %in% c("TopSite")),]$ssFC),color="blue")+
	geom_hline(yintercept=median(allFCmiR1[which(allFCmiR1$site %in% c("NoSite")),]$ssFC),color="red")+
	geom_vline(xintercept=median(allFCmiR1[which(allFCmiR1$site %in% c("Site")),]$halflife_t),color="green")+
	# geom_vline(xintercept=median(allFCmiR1[which(allFCmiR1$site %in% c("TopSite")),]$halflife_t),color="blue")+
	geom_vline(xintercept=median(allFCmiR1[which(allFCmiR1$site %in% c("NoSite")),]$halflife_t),color="red")

# ggsave(plot=p2,file="/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/figures/version_6/FigS7Xpanel2.pdf",width=3,height=3,useDingbats=FALSE)

allFCmiR1$miR = "miR1"
colnames(allFCmiR1)[2]<-"fcSS"
allFCmiR155$miR = "miR155"
allFC <- rbind(allFCmiR1,allFCmiR155)
# allFC <- allFC[allFC$site %in% c("Site","TopSite"),]
# allFCts <- allFC[allFC$site=="TopSite",]
# allFCts$site = "Site"
# allFC <- rbind(allFC,allFCts)

allFCmed <- data.frame(
	site = c("Site","Site","Site","Site"),
	miR = c("miR1","miR1","miR155","miR155"),
	meds = c(median(allFC[allFC$miR=="miR1" & allFC$site=="Site",]$fcSS),median(allFC[allFC$miR=="miR1" & allFC$site=="TopSite",]$fcSS),median(allFC[allFC$miR=="miR155" & allFC$site=="Site",]$fcSS),median(allFC[allFC$miR=="miR155" & allFC$site=="TopSite",]$fcSS))
	)
breaks = 8^(-2:2)
labels = 8^(-2:2)
allFC[which(allFC$halflife_t<2^-5),]$halflife_t=2^-5

p3 <- ggplot(allFC[allFC$miR=="miR1" & allFC$site=="Site",],aes(x=halflife_t,y=fcSS))+geom_point(size=0.5,alpha=0.5)+facet_wrap(~miR+site,ncol=4)+
	scale_x_continuous(trans="log2",expand=c(0,0),limits=2^c(-5,7),breaks=breaks,labels=labels)+theme_tim()+
	scale_y_continuous(limits=c(-2.0,2.0),name=NULL,expand=c(0,0))+
	geom_hline(yintercept=0)+
	geom_hline(data=allFCmed[allFCmed$miR=="miR1",],aes(yintercept=meds),linetype="dashed",color="grey")

p4 <- ggplot(allFC[allFC$miR=="miR155" & allFC$site=="Site",],aes(x=halflife_t,y=fcSS))+geom_point(size=0.5,alpha=0.5)+facet_wrap(~miR+site,ncol=4)+
	scale_x_continuous(trans="log2",expand=c(0,0),limits=2^c(-5,7),breaks=breaks,labels=labels)+theme_tim()+
	scale_y_continuous(limits=c(-2.0,2.0),name=NULL,expand=c(0,0))+
	geom_hline(yintercept=0)+
	geom_hline(data=allFCmed[allFCmed$miR=="miR155",],aes(yintercept=meds),linetype="dashed",color="grey")

# sitemiR1[which(sitemiR1$halflife_t<2^-5),]$halflife_t=2^-5
sitemiR1[which(sitemiR1$TotalCppScore< -2),]$TotalCppScore=-2
sitemiR155[which(sitemiR155$halflife_t<2^-5),]$halflife_t=2^-5

p5 <- ggplot(sitemiR1,aes(x=halflife_t,y=TotalCppScore))+geom_point(size=0.5,alpha=0.5)+
	geom_hline(yintercept=0)+
	scale_y_continuous(limits=c(-2.0,0),name=NULL,expand=c(0,0))+
	scale_x_continuous(trans="log2",expand=c(0,0),limits=2^c(-5,7),breaks=breaks,labels=labels)+
	theme_tim()

p6 <- ggplot(sitemiR155,aes(x=halflife_t,y=TotalCppScore))+geom_point(size=0.5,alpha=0.5)+
	geom_hline(yintercept=0)+
	scale_y_continuous(limits=c(-2.0,0),name=NULL,expand=c(0,0))+
	scale_x_continuous(trans="log2",expand=c(0,0),limits=2^c(-5,7),breaks=breaks,labels=labels)+
	theme_tim()

print(nrow(sitemiR1))
print(nrow(sitemiR155))
print(cor(sitemiR1$halflife_t,sitemiR1$TotalCppScore,method="spearman"))
print(cor(sitemiR155$halflife_t,sitemiR155$TotalCppScore,method="spearman"))
print(cor.test(sitemiR1$halflife_t,sitemiR1$TotalCppScore,method="spearman")[[3]])
print(cor.test(sitemiR155$halflife_t,sitemiR155$TotalCppScore,method="spearman")[[3]])

ggsave(plot=p5,file="/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/figures/version_7/FigS7Xpanel5.pdf",width=1.912,height=1.7)
ggsave(plot=p6,file="/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/figures/version_7/FigS7Xpanel6.pdf",width=1.912,height=1.7)



print("##########################################")
miR1Site = allFC[allFC$miR=="miR1" & allFC$site=="Site",]
# miR1TopSite = allFC[allFC$miR=="miR1" & allFC$site=="TopSite",]
miR155Site = allFC[allFC$miR=="miR155" & allFC$site=="Site",]
# miR155TopSite = allFC[allFC$miR=="miR155" & allFC$site=="TopSite",]
break
print(nrow(miR1Site))
print(nrow(miR1TopSite))
print(nrow(miR155Site))
print(nrow(miR155TopSite))
print(cor(miR1Site$halflife_t,miR1Site$fcSS,method="spearman"))
print(cor(miR1TopSite$halflife_t,miR1TopSite$fcSS,method="spearman"))
print(cor(miR155Site$halflife_t,miR155Site$fcSS,method="spearman"))
print(cor(miR155TopSite$halflife_t,miR155TopSite$fcSS,method="spearman"))
print(cor.test(miR1Site$halflife_t,miR1Site$fcSS,method="spearman")[[3]])
print(cor.test(miR1TopSite$halflife_t,miR1TopSite$fcSS,method="spearman")[[3]])
print(cor.test(miR155Site$halflife_t,miR155Site$fcSS,method="spearman")[[3]])
print(cor.test(miR155TopSite$halflife_t,miR155TopSite$fcSS,method="spearman")[[3]])


ggsave(plot=p3,file="/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/figures/version_7/FigS7Xpanel3.pdf",width=3.4,height=2)
ggsave(plot=p4,file="/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/figures/version_7/FigS7Xpanel4.pdf",width=3.4,height=2)

