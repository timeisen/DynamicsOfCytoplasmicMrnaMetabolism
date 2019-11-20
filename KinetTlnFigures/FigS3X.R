library(ggplot2)
library(plyr)
library(reshape2)
options("scipen"=5)
source("https://raw.githubusercontent.com/briatte/ggcorr/master/ggcorr.R")
source("/lab/solexa_bartel/teisen/RNAseq/Scripts/general/ggplot_theme.R")


stds <- read.table("/lab/solexa_bartel/teisen/RNAseq/Annotation_files/TAIL_seq/list_of_standards_SWE_notation.txt")$V1

m1_TAIL <- read.table("processed_files/miR-1_minus_sample_mean_tails_50tags_background_subtracted_v6.txt")
m155_TAIL <- read.table("processed_files/miR-155_minus_sample_mean_tails_50tags_background_subtracted_v6.txt")
colnames(m1_TAIL) <- c("accession","tp40m.1","tp1hr.1","tp2hr.1","tp4hr.1","tpSS.1")
colnames(m155_TAIL) <- c("accession","tp40m.155","tp1hr.155","tp2hr.155","tp4hr.155","tp8hr.155","tpSS.155")

m1_RPF <- read.table("/lab/solexa_bartel/teisen/RNAseq/kinetics_of_translation/final_analyses/analysis/processed_files/miR-1_RPF_minus_tags_10rpm_cutoff_with_stds_no8hr.txt",head=TRUE)
m155_RPF <- read.table("/lab/solexa_bartel/teisen/RNAseq/kinetics_of_translation/final_analyses/analysis/processed_files/miR-155_RPF_minus_tags_10rpm_cutoff_with_stds.txt",head=TRUE)

colnames(m1_RPF) <- c("accession","tp40m.1","tp1hr.1","tp2hr.1","tp4hr.1","tpSS.1")
colnames(m155_RPF) <- c("accession","tp40m.155","tp1hr.155","tp2hr.155","tp4hr.155","tp8hr.155","tpSS.155")

m1_RNA <- read.table("/lab/solexa_bartel/teisen/RNAseq/kinetics_of_translation/final_analyses/analysis/processed_files/miR-1_RNA_minus_tags_10rpm_cutoff_with_stds_no8hr.txt",head=TRUE)
m155_RNA <- read.table("/lab/solexa_bartel/teisen/RNAseq/kinetics_of_translation/final_analyses/analysis/processed_files/miR-155_RNA_minus_tags_10rpm_cutoff_with_stds.txt",head=TRUE)

colnames(m1_RNA) <- c("accession","tp40m.1","tp1hr.1","tp2hr.1","tp4hr.1","tpSS.1")
colnames(m155_RNA) <- c("accession","tp40m.155","tp1hr.155","tp2hr.155","tp4hr.155","tp8hr.155","tpSS.155")


m_TAIL <- merge(m1_TAIL,m155_TAIL)
m_RPF  <- merge(m1_RPF,m155_RPF)
m_RNA  <- merge(m1_RNA,m155_RNA)

#remove standards
m_TAIL <- m_TAIL[!m_TAIL$accession %in% stds,]
m_RPF <- m_RPF[!m_RPF$accession %in% stds,]
m_RNA <- m_RNA[!m_RNA$accession %in% stds,]

print(nrow(m_TAIL))
print(nrow(m_RPF))
print(nrow(m_RNA))

p1 <- ggcorr(data = NULL, cor_matrix = cor(m_TAIL[, -1], method="spearman"),label=TRUE,label_round = 2,geom="text",hjust=0.75,layout.exp = 1,size = 2,label_size = 2)+theme_tim()

p2 <- ggcorr(data = NULL, cor_matrix = cor(m_RPF[, -1], method="spearman"),label=TRUE,label_round = 2,geom="text",hjust=0.75,layout.exp = 1,size = 2,label_size = 2)+theme_tim()

p3 <- ggcorr(data = NULL, cor_matrix = cor(m_RNA[, -1], method="spearman"),label=TRUE,label_round = 2,geom="text",hjust=0.75,layout.exp = 1,size = 2,label_size = 2)+theme_tim()
break
ggsave(p1,file="figures/version_7/S3X_Test_Tails_cor.pdf",width=4,height=4)
ggsave(p2,file="figures/version_7/S3X_Test_RPF_cor.pdf",width=4,height=4)
ggsave(p3,file="figures/version_7/S3X_Test_RNA_cor.pdf",width=4,height=4)