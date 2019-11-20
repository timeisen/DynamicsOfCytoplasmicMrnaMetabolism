#Figure 3, the heatmaps of individual genes, mRNA paper.
library(ggplot2)
library(plyr)
library(reshape2)
library(scales)
library(cowplot)
args<-commandArgs(trailingOnly=TRUE)
source("ggplot_theme.R")
# theme_tim <- function(base_size = 6, base_family = "Helvetica")
# {
# theme_bw(base_size = base_size, base_family = base_family) +
# theme(
# 	panel.border = element_rect(fill=NA, colour = "black", size=0.5),
# 	panel.grid=element_blank(),
# 	axis.text.y = element_text (size = 6,color="black"),
# 	axis.text.x = element_text (size = 6,color="black"),
# 	# axis.ticks.length = unit(0.25 , "cm"),
# 	#axis.line=element_line(size=.5,color="black"),
# 	#axis.text.x=element_blank(),
#    	#axis.ticks=element_blank(),
#     #axis.title.x=element_blank(),
#     plot.margin=unit(c(0,3,0,0),"mm"),
#     legend.position="none")
# }

bins<-function(df,nt_bins){
	df<-df[order(df[,2]),]
	x <- rep(0:49,each=5)	
	df$x<-x
	x<-unique(x)
	df2<-data.frame(sapply(1:length(x),function(n){return(sum(df[which(df$x==x[n]),][,1]))}),x)
	return(df2)}
bw=5

tp40m<-read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/single_tag_files/background_subtracted_single_tag_files/tp40m_bs_minus_miR-155_v7.txt",sep="\t")
tp1hr<-read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/single_tag_files/background_subtracted_single_tag_files/tp1hr_bs_minus_miR-155_v7.txt",sep="\t")
tp2hr<-read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/single_tag_files/background_subtracted_single_tag_files/tp2hr_bs_minus_miR-155_v7.txt",sep="\t")
tp4hr<-read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/single_tag_files/background_subtracted_single_tag_files/tp4hr_bs_minus_miR-155_v7.txt",sep="\t")
tp8hr<-read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/single_tag_files/background_subtracted_single_tag_files/tp8hr_bs_minus_miR-155_v7.txt",sep="\t")
# tpSS_<-read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/single_tag_files/background_subtracted_single_tag_files/SS_minus_norm_miR-155_v7.txt",sep="\t")
tpSS_<-read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/single_tag_files/background_subtracted_single_tag_files/SS_minus_norm_miR-155_v7_HYBRID20190731.txt",sep = "\t")
# norm_vals <- c(
# 	sum(tp40m[,-252]),
# 	sum(tp1hr[,-252]),
# 	sum(tp2hr[,-252]),
# 	sum(tp4hr[,-252]),
# 	sum(tp8hr[,-252]),
# 	sum(tpSS_[,-252]))

# expr_hl <-read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/processed_files/halflives/miR-155_minus_halflives.txt",head=TRUE)
# expr_hl <- read.table("/lab/solexa_bartel/eichhorn/5EU_RPF/miR-1_miR-155_analysis/Halflife_analysis/Halflife_measurements/miR-155_RNA_minusdox_halflife_SS.txt",sep="\t",head=TRUE)
expr_hl <- read.table("/lab/solexa_bartel/eichhorn/5EU_RPF/miR-1_miR-155_analysis/Halflife_analysis/Halflife_measurements/miR-155_RNA_minusdox_halflife_SS.txt",sep="\t",head=TRUE)
expr_data <- read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/model_input_files/miR-155_minus_sample_mean_tails_50tags_with_expr.txt")
# expr_hl <- expr_hl[,c(1,8,9)]
# expr_hl[,3] <- log(2)/expr_hl[,3]

colnames(expr_hl)[3]<-"hl"
# expr_data <- read.table("/lab/solexa_bartel/eichhorn/5EU_RPF/miR-1_miR-155_analysis/Halflife_analysis/Halflife_measurements/miR-155_RNA_minusdox_halflife_SS.txt",head = TRUE, sep = "\t")
# expr_data <- expr_data[,c(1,8:13)]
# break
# expr_data <- read.table("processed_files/rnaseqDatasets/miR-155_minu_polyA_30minHLnorm.txt",head = TRUE,sep = "\t")
# expr_norm <- colSums(expr_data[,-1])

colnames(expr_data) <- c("accession","tp40m","tp1hr","tp2hr","tp4hr","tp8hr","tpSS_")
colnames(tp40m)[252]<- "accession"
colnames(tp1hr)[252]<- "accession"
colnames(tp2hr)[252]<- "accession"
colnames(tp4hr)[252]<- "accession"
colnames(tp8hr)[252]<- "accession"
colnames(tpSS_)[252]<- "accession"

all_genes<-tp40m[which(
tp40m$accession %in% tp40m$accession &
tp40m$accession %in% tp1hr$accession &
tp40m$accession %in% tp2hr$accession &
tp40m$accession %in% tp4hr$accession &
tp40m$accession %in% tp8hr$accession &
tp40m$accession %in% tpSS_$accession &
tp40m$accession %in% expr_data$accession),]$accession

print(length(all_genes))
expr_hl<-expr_hl[which(expr_hl$accession %in% all_genes),]
annot<-read.table("/lab/solexa_bartel/teisen/RNAseq/Annotation_files/TAIL_seq/symbol_to_accession.txt",head=TRUE)
expr_hl<-merge(expr_hl,annot,by="accession")

# genesets <- rev(c("NM_011580"))
genesets <- rev(c("NM_144797","NM_026032","NM_007907","NM_007393","NM_011580"))

# genesets <- c("NM_001081302",
# 	"NM_001195130",
#    	"NM_007778",
#    	"NM_008086",
#    	"NM_008197",
#    	"NM_008668",
#    	"NM_018871",
#    	"NM_139307",
#    	"NM_009969")

for (i in 1:length(genesets)){
	gene <- genesets[i]
	print(gene)
	print("Halflife:")
	print(expr_hl[which(expr_hl$accession==gene),])
	tp40mg<-tp40m[which(tp40m$accession==gene),]
	tp1hrg<-tp1hr[which(tp1hr$accession==gene),]
	tp2hrg<-tp2hr[which(tp2hr$accession==gene),]
	tp4hrg<-tp4hr[which(tp4hr$accession==gene),]
	tp8hrg<-tp8hr[which(tp8hr$accession==gene),]
	tpSS_g<-tpSS_[which(tpSS_$accession==gene),]
	expr_datag<-expr_data[which(expr_data$accession==gene),]
	norm_vals<-as.numeric(expr_datag[1,-1])

	# norm_vals <- c(44486,112016.4486,161431.2957,288238.6775,598558.9972,5385741.303)


	totals<-nrow(tp40mg)+
	nrow(tp1hrg)+
	nrow(tp2hrg)+
	nrow(tp4hrg)+
	nrow(tp8hrg)+
	nrow(tpSS_g)+
	nrow(expr_datag)
	
	if(totals!=7){break}
	
	print("Found gene in all datasets.")
	
	tp40mg<-data.frame(t(tp40mg[,-252]))
	tp1hrg<-data.frame(t(tp1hrg[,-252]))
	tp2hrg<-data.frame(t(tp2hrg[,-252]))
	tp4hrg<-data.frame(t(tp4hrg[,-252]))
	tp8hrg<-data.frame(t(tp8hrg[,-252]))
	tpSS_g<-data.frame(t(tpSS_g[,-252]))

	# norm_vals <- c(sum(tp40mg[,1]),
	# 	sum(tp1hrg[,1]),
	# 	sum(tp2hrg[,1]),
	# 	sum(tp4hrg[,1]),
	# 	sum(tp8hrg[,1]),
	# 	sum(tpSS_g[,1]))

	# norm_vals <- c(
	# 	sum(tp40mg[-c(1:8,251),1]),
	# 	sum(tp1hrg[-c(1:8,251),1]),
	# 	sum(tp2hrg[-c(1:8,251),1]),
	# 	sum(tp4hrg[-c(1:8,251),1]),
	# 	sum(tp8hrg[-c(1:8,251),1]),
	# 	sum(tpSS_g[-c(1:8,251),1]))

	norm_vals <- c(
		sum(tp40mg[,1]),
		sum(tp1hrg[,1]),
		sum(tp2hrg[,1]),
		sum(tp4hrg[,1]),
		sum(tp8hrg[,1]),
		sum(tpSS_g[,1]))

	print(norm_vals)
	
	# norm_vals <- rev(c(1,
	# 1.74449207,
	# 0.193279838,
	# 0.12829349,
	# 0.079662677,
	# 0.048435329))

	tp40mg[,1]<-tp40mg[,1]/sum(tp40mg[,1])
	tp1hrg[,1]<-tp1hrg[,1]/sum(tp1hrg[,1])
	tp2hrg[,1]<-tp2hrg[,1]/sum(tp2hrg[,1])
	tp4hrg[,1]<-tp4hrg[,1]/sum(tp4hrg[,1])
	tp8hrg[,1]<-tp8hrg[,1]/sum(tp8hrg[,1])
	tpSS_g[,1]<-tpSS_g[,1]/sum(tpSS_g[,1])

		
	tp40mg$tail_length<-0:250
	tp1hrg$tail_length<-0:250
	tp2hrg$tail_length<-0:250
	tp4hrg$tail_length<-0:250
	tp8hrg$tail_length<-0:250
	tpSS_g$tail_length<-0:250
	
	tp40mg <- tp40mg[-251,]
	tp1hrg <- tp1hrg[-251,]
	tp2hrg <- tp2hrg[-251,]
	tp4hrg <- tp4hrg[-251,]
	tp8hrg <- tp8hrg[-251,]
	tpSS_g <- tpSS_g[-251,]
	
	tp40mg<-bins(tp40mg,bw)
	tp1hrg<-bins(tp1hrg,bw)
	tp2hrg<-bins(tp2hrg,bw)
	tp4hrg<-bins(tp4hrg,bw)
	tp8hrg<-bins(tp8hrg,bw)
	tpSS_g<-bins(tpSS_g,bw)

	#these next two blocks are to create the correct horizontal lines at the end of geom_step plotting
	tp40mg <- rbind(tp40mg,tp40mg[50,])
	tp1hrg <- rbind(tp1hrg,tp1hrg[50,])
	tp2hrg <- rbind(tp2hrg,tp2hrg[50,])
	tp4hrg <- rbind(tp4hrg,tp4hrg[50,])
	tp8hrg <- rbind(tp8hrg,tp8hrg[50,])
	tpSS_g <- rbind(tpSS_g,tpSS_g[50,])

	tp40mg[50,2]<-50
	tp1hrg[50,2]<-50
	tp2hrg[50,2]<-50
	tp4hrg[50,2]<-50
	tp8hrg[50,2]<-50
	tpSS_g[50,2]<-50

	print("Finished binning ...")
	colnames(tp40mg)<-c("count","tail_length")
	colnames(tp1hrg)<-c("count","tail_length")
	colnames(tp2hrg)<-c("count","tail_length")
	colnames(tp4hrg)<-c("count","tail_length")
	colnames(tp8hrg)<-c("count","tail_length")
	colnames(tpSS_g)<-c("count","tail_length")
	
	tp40mg$sample<-"40 min"
	tp1hrg$sample<-"1 h"
	tp2hrg$sample<-"2 h"
	tp4hrg$sample<-"4 h"
	tp8hrg$sample<-"8 h"
	tpSS_g$sample<-"SS"
	
	print(norm_vals)
	tp40mg$norm<-1/norm_vals[1]
	tp1hrg$norm<-1/norm_vals[2]
	tp2hrg$norm<-1/norm_vals[3]
	tp4hrg$norm<-1/norm_vals[4]
	tp8hrg$norm<-1/norm_vals[5]
	tpSS_g$norm<-1/norm_vals[6]
	#tp40m$count[tp40m$count<0] = 0
	#tp1hr$count[tp1hr$count<0] = 0
	#tp2hr$count[tp2hr$count<0] = 0
	#tp4hr$count[tp4hr$count<0] = 0
	#tp8hr$count[tp8hr$count<0] = 0
	#tpSS_$count[tpSS_$count<0] = 0
	
	tp40mg$norm1<-sum(tp40mg$count)
	tp1hrg$norm1<-sum(tp1hrg$count)
	tp2hrg$norm1<-sum(tp2hrg$count)
	tp4hrg$norm1<-sum(tp4hrg$count)
	tp8hrg$norm1<-sum(tp8hrg$count)
	tpSS_g$norm1<-sum(tpSS_g$count)
	
	
	# tp40min$counts<-809790
	# tp1hr$counts<-1412688
	# tp2hr$counts<-313036
	# tp4hr$counts<-259737
	# tp8hr$counts<-322556
	# tpSS$counts<-140468
	
	
	
	cols<-c(
	"black",
	"#8856a7",
	"#2b8cbe",
	"#31a354",
	"#fec44f",
	"#f03b20")
	
	tails<-NULL
	
	
	tails<-rbind(tp40mg,tp1hrg,tp2hrg,tp4hrg,tp8hrg,tpSS_g)
	#tails<-rbind(tp40min,tp1hr)
	# tails[which(tails$sample=="SS" & tails$tail_length==49),]$count = max(tails$count)
	# tails[which(tails$sample=="SS" & tails$tail_length==39),]$count = max(tails$count)
	# tails[which(tails$sample=="SS" & tails$tail_length==29),]$count = max(tails$count)
	# tails[which(tails$sample=="SS" & tails$tail_length==19),]$count = max(tails$count)
	# tails[which(tails$sample=="SS" & tails$tail_length==9),]$count = max(tails$count)
	tails$count = tails$count/sum(tails$count/tails$norm)


	tails$sample<-factor(tails$sample,levels=c("SS","8 h","4 h","2 h","1 h","40 min"))
	tails$tail_length<-as.numeric(as.character(tails$tail_length))
	tails$rescale<-0
	tails$rescale <- as.numeric(tails$count)/sum(as.numeric(tails$count))*10
	
	
	print(median(tails$count))
	print(head(tails))
	print(min(tails$count))
	tails$offset<-0
	
	tails[which(tails$sample=="SS"),]$offset <- tails[which(tails$sample=="SS"),]$rescale
	tails[which(tails$sample=="8 h"),]$offset <- tails[which(tails$sample=="8 h"),]$rescale+500
	tails[which(tails$sample=="4 h"),]$offset <- tails[which(tails$sample=="4 h"),]$rescale+400
	tails[which(tails$sample=="2 h"),]$offset <- tails[which(tails$sample=="2 h"),]$rescale+300
	tails[which(tails$sample=="1 h"),]$offset <- tails[which(tails$sample=="1 h"),]$rescale+200
	tails[which(tails$sample=="40 min"),]$offset <- tails[which(tails$sample=="40 min"),]$rescale+100
	
	#tails$offset <-tails$rescale+100*(as.numeric(tails$sample)-1)
	scalerange <- as.numeric(quantile(tails$rescale,p=c(0.05,.95)))
	
	print(scalerange)
	gradientends <- scalerange + rep(c(0,5,4,3,2,1)*100, each=2)
	
	cols_hm<-cols
	colorends <- c("white", cols_hm[1], "white", cols_hm[2], "white", cols_hm[3],"white",cols_hm[4],"white",cols_hm[5],"white",cols_hm[6])
	print("Plotting ...")


	p1 <- ggplot(tails, aes(x=tail_length*5+2.5, y=sample))+geom_tile(aes(fill = offset))+
	scale_fill_gradientn(colours = colorends, values = rescale(gradientends))+
	scale_x_continuous(name=NULL,limits=c(10,250),breaks=c(10,50,100,150,200,250)) + 
	scale_y_discrete(name=NULL,expand=c(0,0))+theme_tim()+
	theme(
		legend.position="none",
		axis.ticks.y=element_blank(),
		axis.line=element_blank(),
		axis.text.y=element_blank()) + 
	geom_segment(aes(x = 10, xend = 250, y = 0, yend = 0))

	p2<-ggplot(tails,aes(x=tail_length*5,y=count/norm,color=sample))+geom_step(size=0.5,direction="hv")+scale_colour_manual(values=cols)+
	scale_y_continuous(name=NULL,expand = c(0,0),limits=c(NA,max(tails$count/tails$norm+tails$count/tails$norm/10)))+theme_tim()+
	scale_x_continuous(name=NULL,limits=c(10,250))+
	theme(
		legend.position="none",
		axis.text.x=element_blank(),
		axis.line.x=element_blank(),
	   	axis.ticks.x=element_blank())+
	geom_segment(aes(x = 10, xend = 250, y = 0, yend = 0),size = 0.2)

	
	
	# p3<-ggdraw() +
	#   draw_plot(p2, 0, .25, 1, .75) +
	#   draw_plot(p1, 0, 0, 1, .25)
	
	p3 <- plot_grid(p2,p1,align='v',labels=NULL,nrow=2,rel_heights=c(3,1))
	ggsave(plot=p3,file=paste0("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/figures/version_9/Fig3Cpanel",gene,".pdf"),width=2.2,height=2.2)
}