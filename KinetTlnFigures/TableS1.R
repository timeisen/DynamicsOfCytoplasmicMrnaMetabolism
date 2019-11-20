#This script generates the first table.
#Format should be: alpha miR-1 uninduced, alpha, st, dead; miR-155 uninduced alpha, st, dead, lifetime,
library(plyr)
library(tidyverse)

pLogEl155 = plogis(230, loc = 263.95156,scale = 11.05133)
pLogEl1 = plogis(230, loc = 265.44891,scale = 13.97311)

miR155_HL <- read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/processed_files/halflifeComparisons/miR-155_minu_PAL_halflives_logspace_global_offset_SS_RNAseqScaling.txt",head=TRUE)
miR155_HL <- miR155_HL[,c(1,5)]
colnames(miR155_HL)[2] <- "miR155_HL"
miR155_HL[miR155_HL$miR155_HL < 0.1,]$miR155_HL = 0.1
miR155_HL[miR155_HL$miR155_HL > 1E2,]$miR155_HL = 1E2

paste("miR155 half-lives, number of genes:", nrow(miR155_HL))

miR1_HL <- read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/processed_files/halflifeComparisons/miR-1_minu_PAL_halflives_logspace_global_offset_SS_RNAseqScaling.txt",head=TRUE)
miR1_HL <- miR1_HL[,c(1,5)]
colnames(miR1_HL)[2] <- "miR1_HL"
miR1_HL[miR1_HL$miR1_HL < 0.1,]$miR1_HL = 0.1
miR1_HL[miR1_HL$miR1_HL > 1E2,]$miR1_HL = 1E2

paste("miR1 half-lives, number of genes:",nrow(miR1_HL))

miR155_minus_RATES <- read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/rate_constant_measurements/miR-155_minus_samples_UNLINKV3_V81H5_run3_high_precision_reformat.txt")
colnames(miR155_minus_RATES)<-c("accession","st_m155_minus","a_m155_minus","k_m155_minus","b_m155_minus","r_m155_minus")

miR155_minus_RATES[miR155_minus_RATES$k_m155_minus < 0.03,]$k_m155_minus = 0.03
miR155_minus_RATES[miR155_minus_RATES$k_m155_minus > 30,]$k_m155_minus = 30
miR155_minus_RATES[miR155_minus_RATES$b_m155_minus*pLogEl155 < 0.003,]$b_m155_minus = 0.003/pLogEl155
miR155_minus_RATES[miR155_minus_RATES$b_m155_minus*pLogEl155 > 3,]$b_m155_minus = 3/pLogEl155
miR155_minus_RATES[miR155_minus_RATES$a_m155_minus < 1E-8,]$a_m155_minus = 1E-8
miR155_minus_RATES[miR155_minus_RATES$a_m155_minus > 1E-5,]$a_m155_minus = 1E-5


paste("miR155 model rates, number of genes:",nrow(miR155_minus_RATES))

# miR155_plus_RATES <- read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/rate_constant_measurements/miR-155_plus_samples_v74_plogis_linked_alpha_run3.txt")
# miR155_plus_RATES <- miR155_plus_RATES[,c(1,4,5)]
# colnames(miR155_plus_RATES) <- c("accession","k_m155_plus","r_m155_plus")

miR1_minus_RATES <- read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/rate_constant_measurements/miR-1_minus_samples_UNLINKV3_V81H5_run3_high_precision_reformat.txt")
colnames(miR1_minus_RATES)<-c("accession","st_m1_minus","a_m1_minus","k_m1_minus","b_m1_minus","r_m1_minus")

miR1_minus_RATES[miR1_minus_RATES$k_m1_minus < 0.03,]$k_m1_minus = 0.03
miR1_minus_RATES[miR1_minus_RATES$k_m1_minus > 30,]$k_m1_minus = 30
miR1_minus_RATES[miR1_minus_RATES$b_m1_minus*pLogEl1 < 0.003,]$b_m1_minus = 0.003/pLogEl1
miR1_minus_RATES[miR1_minus_RATES$b_m1_minus*pLogEl1 > 3,]$b_m1_minus = 3/pLogEl1
# miR1_minus_RATES[miR1_minus_RATES$a_m1_minus < 1E-8,]$a_m1_minus = 1E-8
miR1_minus_RATES[miR1_minus_RATES$a_m1_minus > 1E-5,]$a_m1_minus = 1E-5

paste("miR1 model rates, number of genes:",nrow(miR1_minus_RATES))
# miR1_plus_RATES <- read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/rate_constant_measurements/miR-1_plus_samples_v74_plogis_linked_run3.txt")
# miR1_plus_RATES <- miR1_plus_RATES[,c(1,4,5)]
# colnames(miR1_plus_RATES) <- c("accession","k_m1_plus","r_m1_plus")

miR155_nosite <- read.table("/archive/bartel/2015.03.04-11263/solexa_bartel/eichhorn/Ribosome_profiling_paper/Sites/Mouse/miR-155/mmu-miR-155_NoSite.txt")
miR155_site <-read.table("/archive/bartel/2015.03.04-11263/solexa_bartel/eichhorn/Ribosome_profiling_paper/Sites/Mouse/miR-155/mmu-miR-155_tpUTR_oneormore.txt")
miR155_topsite <-read.table("/lab/solexa_bartel/eichhorn/5EU_paper/Eichhorn_data/RESTORE/miR-155_3T3_repressed_genes_25pct_down.txt")

miR1_nosite <- read.table("/archive/bartel/2015.03.04-11263/solexa_bartel/eichhorn/Ribosome_profiling_paper/Sites/Mouse/miR-1/mmu-miR-1_NoSite.txt")
miR1_site <- read.table("/archive/bartel/2015.03.04-11263/solexa_bartel/eichhorn/Ribosome_profiling_paper/Sites/Mouse/miR-1/mmu-miR-1_tpUTR_oneormore.txt")
miR1_topsite <- read.table("/lab/solexa_bartel/eichhorn/5EU_paper/Eichhorn_data/miR-1_3T3_repressed_genes_25pct_down.txt")


ANNOT <- read.table("/lab/solexa_bartel/teisen/RNAseq/Annotation_files/TAIL_seq/symbol_to_accession.txt",head=TRUE)

list_m <- list(
	miR1_minus_RATES,
	# miR1_plus_RATES,
	miR1_HL,
	miR155_minus_RATES,
	# miR155_plus_RATES,
	miR155_HL
	)

tab_m <- join_all(list_m,by="accession",type="full")
tab_m$miR1_site <- NA
tab_m[tab_m$accession %in% miR1_nosite$V1,]$miR1_site <- "No Site, miR-1"
tab_m[tab_m$accession %in% miR1_site$V1,]$miR1_site <- "Site, miR-1"
tab_m[tab_m$accession %in% miR1_topsite$V1,]$miR1_site <- "Top Site, miR-1"
tab_m$miR155_site <- NA
tab_m[tab_m$accession %in% miR155_nosite$V1,]$miR155_site <- "No Site, miR-155"
tab_m[tab_m$accession %in% miR155_site$V1,]$miR155_site <- "Site, miR-155"
tab_m[tab_m$accession %in% miR155_topsite$V1,]$miR155_site <- "Top Site, miR-155"

tab_m <- merge(tab_m,ANNOT,by="accession")

tab_m <- tab_m[,c(1,16,2:5,7,14,8:11,13,15)]
print(cor(tab_m[,c(3:7,9:13)],method = 'spearman', use = 'pairwise.complete.obs'))

## Cell line 1: miR-155
## Cell line 2: miR-1
tab_m <- rbind(c(rep(NA,2),rep("Cell Line 2 (miR-1)",6),rep("Cell Line 1 (miR-155)",6)),tab_m)

colnames(tab_m) <- c(
	"AccessionCode",
	"GeneSymbol",
	"StartingTailLength",
	"ProductionRate",
	"DeadenylationRate",
	"DecappingRate",
	"HalfLife",
	"SiteType",
	"StartingTailLength",
	"ProductionRate",
	"DeadenylationRate",
	"DecappingRate",
	"HalfLife",
	"SiteType"
	)

write.table(tab_m, "/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/figures/version_9/TableS1.txt",
            na = "--",
            row.names = FALSE,
            col.names = TRUE,
            append = FALSE,
            sep = "\t",
            quote=FALSE)

