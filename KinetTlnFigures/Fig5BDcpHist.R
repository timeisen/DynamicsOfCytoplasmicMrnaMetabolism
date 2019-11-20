library(data.table)
library(spatstat)
library(ggplot2)
source("/lab/solexa_bartel/teisen/RNAseq/Scripts/general/ggplot_theme.R")

annot<-read.table("/lab/solexa_bartel/teisen/RNAseq/Annotation_files/TAIL_seq/symbol_to_accession.txt",head=TRUE)
rates_minus<-read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/rate_constant_measurements/miR-155_minus_samples_UNLINKV3_V81H5_run3_high_precision_reformat.txt")
colnames(rates_minus)<-c("accession","st","a","k","b","residual")
hl <- read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/processed_files/halflifeComparisons/miR-155_minu_PAL_halflives_logspace_global_offset_SS_RNAseqScaling.txt",head=TRUE,sep="\t")
rates_minus<-merge(rates_minus,hl,by="accession")
rates_minus<-rates_minus[order(rates_minus$halflife),]

#Generates Dcp Rates for endogenous genes
all_genes_df = fread('/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/figures/other/all_dcp_rates_v8.txt',head=TRUE,sep="\t")
all_genes_df <- all_genes_df[which(all_genes_df$time==5965),]
all_genes_df$accession <- rates_minus$accession
all_genes_df[,231:250] <- apply(all_genes_df[,231:250],1,mean)
all_genes_df$mean <- apply(all_genes_df[,1:250],1,function(x){weighted.mean(249:0,x)})
all_genes_df <- merge(all_genes_df,annot,by="accession")
all_genes_df$sample <- "Data"

dcp_mean <- data.frame(all_genes_df$symbol,all_genes_df$mean)
colnames(dcp_mean) <- c("symbol","dcp_mean_tl")

write.table(dcp_mean,file = "processed_files/DcpMeanTlmiR-155_v7_data.txt",quote = FALSE, sep = "\t",row.names = FALSE)

#Generate Dcp Rates for simulated NMD targets. 
# nmd_sim_targets <- read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/figures/other/all_dcp_rates_NMD_simulation_v1.txt",head=TRUE,sep="\t")
# colnames(nmd_sim_targets) <- c(249:0,"time","gene")
# nmd_sim_targets <- nmd_sim_targets[which(nmd_sim_targets$time == 5965),]
# nmd_sim_targets$accession <- nmd_sim_targets$gene
# nmd_sim_targets$mean <- apply(nmd_sim_targets[,1:250],1,function(x){weighted.mean(249:0,x)})
# nmd_sim_targets$symbol <- nmd_sim_targets$accession
# nmd_sim_targets$sample <- "NMD simulated"

# all_genes_df <- rbind(all_genes_df,nmd_sim_targets)

# upf1 <- read.table("/lab/solexa_bartel/teisen/RNAseq/Annotation_files/nmd_targets_hurt_2013.txt",head=TRUE,sep="\t",skip=1)
# upf1 <- upf1[order(-upf1$Geometric.mean),]
# colnames(upf1)[1] <- "symbol"
all_genes_df <- all_genes_df[order(-all_genes_df$mean),]
all_genes_df$rank = 1:nrow(all_genes_df)


# m <- merge(all_genes_df[,-c(2:251)],upf1,by="symbol",all=TRUE)
# NMDgenes <- read.table("other_analyses/NMDsimulations/NMDgenesLykkeAndsersen2014.txt",head=TRUE,sep="\t")
# NMDgenesList <- NMDgenes[which(NMDgenes$subtype == "decap_NMD"),]$symbol
# all_genes_df$nmd <- "False"
# all_genes_df[which(all_genes_df$symbol %in% NMDgenesList),]$nmd <- "true"

# maxRankSim <- max(all_genes_df[all_genes_df$sample == "NMD simulated",]$rank)
NMDcalls <- all_genes_df[all_genes_df$rank < 49,c(1,253:257)]
print(NMDcalls)

# p1 <- ggplot(all_genes_df,aes(x =rank, y = mean, labels = symbol,color=sample))+
# 	geom_point(size = 0.5)+
# 	geom_text(data=all_genes_df[which(all_genes_df$rank < maxRankSim & all_genes_df$sample == "Data"),],aes(x =rank, y = mean,label=symbol),size = 1)+
# 	theme_tim()
NMDgenes <- 	c("Gadd45g",
				 "Mat2a",
				 "Serpine1",
				 "H2afx",
				 "Gadd45b",
				 "Hnrnph3")
NMDcalls <- NMDcalls[which(NMDcalls$symbol %in% NMDgenes),]



# all_genes_df[all_genes_df$mean < 20,]$mean = 20
p2 <- ggplot(all_genes_df,aes(x = mean, y = ..density..))+
	stat_bin(geom="step", bins = 250,position='identity')+
	scale_x_continuous(expand=c(0,0),limits=c(0,250),breaks=c(0,50,100,150,200,250))+
	scale_y_continuous(expand=c(0,0))+
	# geom_vline(xintercept = NMDcalls$mean)+
	theme_tim()

# all_genes_df <- all_genes_df[order(-all_genes_df$mean),]
# all_genes_df$ecdf <- c(all_genes_df$mean[1],sapply(2:nrow(all_genes_df),function(x){sum(all_genes_df$mean[1:x])}))
# all_genes_df$ecdf <- all_genes_df$ecdf/max(all_genes_df$ecdf)

obj_all_genes_df_ecdf = ecdf(-all_genes_df$mean)
all_genes_df_ecdf <- data.frame(tail_length = 0:250, Valecdf = obj_all_genes_df_ecdf(-c(0:250)))
p2a <- ggplot(all_genes_df_ecdf,aes(x = tail_length, y = Valecdf))+
	geom_step()+
	scale_x_continuous(expand=c(0,0),limits=c(0,250),breaks=c(0,50,100,150,200,250))+
	scale_y_continuous(expand=c(0,0),limits = c(0,1))+
	# geom_vline(xintercept = NMDcalls$mean)+
	theme_tim()

p2inset <- ggplot(all_genes_df,aes(x = mean, y = ..density..))+
	stat_bin(geom="step", bins = 100,position='identity')+
	scale_x_continuous(expand=c(0,0),breaks=c(55,100,120))+
	scale_y_continuous(expand=c(0,0),breaks = c(0,0.001))+
	coord_cartesian(xlim=c(55,120),ylim = c(0,0.001))+
	geom_vline(xintercept = NMDcalls$mean,color="grey",linetype = "dashed")+
	theme_tim()+
	geom_segment(aes(x = 55 , xend = 55, y = 0, yend = 0.001),color="black")+
	theme(axis.line.y = element_blank())


all_genes_quant <- all_genes_df[which(all_genes_df$sample == "Data"),]

print("percentages of mean tl < 50 and < 25, respectively:")
print(1 - obj_all_genes_df_ecdf(-50))
print(1 - obj_all_genes_df_ecdf(-25))

ggsave(plot=p2,file="figures/version_9/Fig5BDcpHist.pdf",useDingbats = FALSE,width = 2,height = 2)
ggsave(plot=p2a,file="figures/version_9/Fig5BDcpHistv2.pdf",useDingbats = FALSE,width = 2,height = 2)
ggsave(plot=p2inset,file="figures/version_9/Fig5BDcpHistInset.pdf",useDingbats = FALSE,width = 1.25,height = 1)
