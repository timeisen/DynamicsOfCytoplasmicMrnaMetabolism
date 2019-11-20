library(tidyverse)
library(data.table)
source("/lab/solexa_bartel/teisen/RNAseq/Scripts/general/ggplot_theme.R")

# bins<-function(df,nt_bins){
# 	df<-df[order(df[,1]),]
# 	x <- rep(0:(250/nt_bins-1),each=nt_bins)	
# 	df$x<-x
# 	x<-unique(x)
# 	df2<-data.frame(sapply(1:length(x),function(n){return(sum(df[which(df$x==x[n]),][,2]))}),x)
# 	return(df2)}

all_genes_df = fread('/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/figures/other/all_dcp_rates_v8.txt',head=TRUE,sep="\t")
# all_genes_df_norm = fread('/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/figures/other/all_dcp_rates_v2.txt',head=TRUE,sep="\t")

# all_genes_df <- all_genes_df[which(all_genes_df$time==5967 & all_genes_df$gene != 7),]
all_genes_df <- all_genes_df[which(all_genes_df$time==5965),]
# all_genes_df[,1:250] <- all_genes_df[,1:250]/rowSums(all_genes_df[,1:250])
all_genes_df_test <- all_genes_df

# print(nrow(all_genes_df[log10(rowSums(all_genes_df[,1:250]))> 0,]))
# all_genes_df <- all_genes_df[log10(rowSums(all_genes_df[,1:250]))< -0,]
# all_genes_df <- all_genes_df[-c(3),]

ngene = 100
all_genes_df_1 <- colSums(all_genes_df[sample(1:nrow(all_genes_df),ngene),])
all_genes_df_2 <- colSums(all_genes_df[sample(1:nrow(all_genes_df),ngene),])
all_genes_df_3 <- colSums(all_genes_df[sample(1:nrow(all_genes_df),ngene),])
all_genes_df_4 <- colSums(all_genes_df[sample(1:nrow(all_genes_df),ngene),])
all_genes_df_5 <- colSums(all_genes_df[sample(1:nrow(all_genes_df),ngene),])
all_genes_df_6 <- colSums(all_genes_df[sample(1:nrow(all_genes_df),ngene),])
all_genes_df_7 <- colSums(all_genes_df[sample(1:nrow(all_genes_df),ngene),])
all_genes_df_8 <- colSums(all_genes_df[sample(1:nrow(all_genes_df),ngene),])
all_genes_df_9 <- colSums(all_genes_df[sample(1:nrow(all_genes_df),ngene),])
all_genes_df_X <- colSums(all_genes_df[,1:250])

# hista = colSums(all_genes_df)
# histanorm = colSums(all_genes_df_norm[which(all_genes_df_norm$time==5967),])
print(length(all_genes_df_X))
hm <- data.frame(
	tl = 249:0,
	all_genes_1 = all_genes_df_1[1:250],
	all_genes_2 = all_genes_df_2[1:250],
	all_genes_3 = all_genes_df_3[1:250],
	all_genes_4 = all_genes_df_4[1:250],
	all_genes_5 = all_genes_df_5[1:250],
	all_genes_6 = all_genes_df_6[1:250],
	all_genes_7 = all_genes_df_7[1:250],
	all_genes_8 = all_genes_df_8[1:250],
	all_genes_9 = all_genes_df_9[1:250])

hm2 <- data.frame(
	tl = 249:0,
	all_genes_X = all_genes_df_X[1:250]/sum(all_genes_df_X[1:250]))
hm3 <- data.frame(
	tl = 249:0,
	all_genes_X = all_genes_df_X[1:250]/sum(all_genes_df_X[1:250]))

hm3$ecdf <- c(hm3$all_genes_X[1],sapply(2:nrow(hm3),function(x){sum(hm3$all_genes_X[1:x])}))

print("decapping at or below 50, percentage:")
print(1 - hm3[which(hm3$tl == 50),]$ecdf)
print("decapping at or below 25, percentage:")
print(1 - hm3[which(hm3$tl == 25),]$ecdf)


hmm <- melt(hm,id="tl")
p3 <- ggplot(hm2[which(hm2$tl!=0),], aes(x = tl, y = all_genes_X)) +geom_step()+scale_x_continuous()+theme_bw() #histogram
p3a <- ggplot(hm2, aes(x = tl, y = all_genes_X)) +geom_step()+scale_x_continuous()+theme_bw() #histogram

p4 <- ggplot(hmm[which(hmm$tl!=0),], aes(x = tl, y = value,color=variable)) +geom_step()+scale_x_continuous()+theme_bw() #histogram

	# all_genes_norm = histanorm[1:250]/sum(histanorm[1:250]))
p1 <- ggplot(hm3, aes(x = tl, y = all_genes_X)) +geom_step()+scale_x_continuous(expand=c(0,0),limits=c(0,250),breaks=c(0,50,100,150,200,250))+scale_y_continuous(name=NULL,expand=c(0,0),limits = c(0,0.03))+theme_tim() #histogram

p1c <- ggplot() +
	geom_step(data = hm3[hm3$tl >= 20,], aes(x = tl, y = all_genes_X)) +
	geom_step(data = hm3[hm3$tl < 20,], aes(x = tl, y = all_genes_X), linetype = 'dashed') +
	scale_x_continuous(expand=c(0,0),limits=c(0,250),breaks=c(0,50,100,150,200,250)) +
	scale_y_continuous(name=NULL,expand=c(0,0),limits = c(0,0.03)) +
	theme_tim() #histogram


print("dataset max:")
print(max(hm3$all_genes_X))

p2 <- ggplot(hm3, aes(x = tl, y = ecdf)) +geom_line()+scale_x_continuous(expand=c(0,0),limits=c(0,250))+scale_y_continuous(name=NULL,expand=c(0,0))+theme_tim() #ecdf

hmLast20 <- data.frame(t(replicate(20,apply(hm3[231:250,],2,mean))))
hmLast20$tl <- 19:0
hm4 <- rbind(hm3[1:230,],hmLast20)
hm4$ecdf <- c(hm4$all_genes_X[1],sapply(2:nrow(hm4),function(x){sum(hm4$all_genes_X[1:x])}))


p1a <- ggplot(hm4, aes(x = tl, y = all_genes_X)) +geom_step()+scale_x_continuous(expand=c(0,0),limits=c(0,250),breaks=c(0,50,100,150,200,250))+scale_y_continuous(name=NULL,expand=c(0,0))+theme_tim() #histogram

p1b <- ggplot(hm3, aes(x = tl, y = ecdf)) +geom_step()+scale_x_continuous(expand=c(0,0),limits=c(0,250),breaks=c(0,50,100,150,200,250))+scale_y_continuous(name=NULL,expand=c(0,0),limits = c(0,1))+theme_tim() #histogram

# hmmelt <- melt(hm,id="tl")
# p2 <- ggplot(hmmelt, aes(x = tl, y = value,color=variable)) +geom_step()+scale_x_continuous(expand=c(0,0),limits=c(0,250))+scale_y_continuous(name=NULL,expand=c(0,0))+theme_tim()+theme(legend.position="right") #histogram

# hm$ecdf <- c(hm$all_genes[1],sapply(1:nrow(hm),function(x){sum(hm$all_genes[1:x])}))
# ggsave(plot=p2,file="/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/figures/other/DcpPositionmRNAsvsGenes.pdf",width=5,height=5)
print(hm[which.max(hm$all_genes),])

ggsave(plot=p1,file="/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/figures/version_9/Fig4E.pdf",width=2,height=2)
ggsave(plot=p1c,file="/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/figures/version_9/Fig4Ev4.pdf",width=2,height=2)
ggsave(plot=p1a,file="/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/figures/version_9/Fig4Ev2.pdf",width=2,height=2)
ggsave(plot=p1b,file="/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/figures/version_9/Fig4Ev3.pdf",width=2,height=2)