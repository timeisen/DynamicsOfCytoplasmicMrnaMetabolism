#Fig7BC mRNA paper. 
#Single tails and normalized single tails
library(plyr)
library(tidyverse)
library(spatstat)
source("ggplot_theme.R")
cols<-rev(c(
"#8856a7",
"#2b8cbe",
"#31a354",
"#fec44f",
"#f03b20"))

bins<-function(vec,nt_bins){
	vec <- lapply(0:(121),function(x){
		if(x == 121){
			return(vec[243])
		} else {
			return(mean(vec[(nt_bins*x+1):((x+1)*nt_bins)]))
		}
		})
	return(as.numeric(vec))}

bw=2

stds<-read.table("/lab/solexa_bartel/teisen/RNAseq/Annotation_files/TAIL_seq/list_of_standards_SWE_notation.txt")
print("Reading data files ...")
actD_0hr <- read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-1_actD_samples/analysis_20171030/model_input_files/global_normalized/actD_0hr_minus_background_subtracted_single_tag_v6.txt",sep="\t")
actD_1hr <- read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-1_actD_samples/analysis_20171030/model_input_files/global_normalized/actD_1hr_minus_background_subtracted_single_tag_v6.txt",sep="\t")
actD_3hr <- read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-1_actD_samples/analysis_20171030/model_input_files/global_normalized/actD_3hr_minus_background_subtracted_single_tag_v6.txt",sep="\t")
actD_7hr <- read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-1_actD_samples/analysis_20171030/model_input_files/global_normalized/actD_7hr_minus_background_subtracted_single_tag_v6.txt",sep="\t")
actD_15h <- read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-1_actD_samples/analysis_20171030/model_input_files/global_normalized/actD_15h_minus_background_subtracted_single_tag_v6.txt",sep="\t")

actD_0hr_m <- actD_0hr[-252]/sum(actD_0hr[-252])
actD_1hr_m <- actD_1hr[-252]/sum(actD_1hr[-252])
actD_3hr_m <- actD_3hr[-252]/sum(actD_3hr[-252])
actD_7hr_m <- actD_7hr[-252]/sum(actD_7hr[-252])
actD_15h_m <- actD_15h[-252]/sum(actD_15h[-252])

actD_0hr_m[252] <- actD_0hr[252] 
actD_1hr_m[252] <- actD_1hr[252] 
actD_3hr_m[252] <- actD_3hr[252] 
actD_7hr_m[252] <- actD_7hr[252] 
actD_15h_m[252] <- actD_15h[252] 

actD_0hr <- actD_0hr_m
actD_1hr <- actD_1hr_m
actD_3hr <- actD_3hr_m
actD_7hr <- actD_7hr_m
actD_15h <- actD_15h_m

print("Gathering data ...")
colnames(actD_0hr)[252]<-"accession"
colnames(actD_1hr)[252]<-"accession"
colnames(actD_3hr)[252]<-"accession"
colnames(actD_7hr)[252]<-"accession"
colnames(actD_15h)[252]<-"accession"

actD_0hr <- actD_0hr[which(!actD_0hr$accession %in% stds$V1),]
actD_1hr <- actD_1hr[which(!actD_1hr$accession %in% stds$V1),]
actD_3hr <- actD_3hr[which(!actD_3hr$accession %in% stds$V1),]
actD_7hr <- actD_7hr[which(!actD_7hr$accession %in% stds$V1),]
actD_15h <- actD_15h[which(!actD_15h$accession %in% stds$V1),]

m <- data.frame(
	tp0 = bins(colSums(actD_0hr[,-c(1:8,252)]),bw),
	tp1 = bins(colSums(actD_1hr[,-c(1:8,252)]),bw),
	tp3 = bins(colSums(actD_3hr[,-c(1:8,252)]),bw),
	tp7 = bins(colSums(actD_7hr[,-c(1:8,252)]),bw),
	t15 = bins(colSums(actD_15h[,-c(1:8,252)]),bw))

print("Calculating mean tails ...")
# actD_0hr_t <- actD_0hr[,-252]
# actD_1hr_t <- actD_1hr[,-252]
# actD_3hr_t <- actD_3hr[,-252]
# actD_7hr_t <- actD_7hr[,-252]
# actD_15h_t <- actD_15h[,-252]

# actD_0hr_t[actD_0hr_t < 0] = 0
# actD_1hr_t[actD_1hr_t < 0] = 0
# actD_3hr_t[actD_3hr_t < 0] = 0
# actD_7hr_t[actD_7hr_t < 0] = 0
# actD_15h_t[actD_15h_t < 0] = 0

# actD_0hr[,-252] <- actD_0hr_t
# actD_1hr[,-252] <- actD_1hr_t
# actD_3hr[,-252] <- actD_3hr_t
# actD_7hr[,-252] <- actD_7hr_t
# actD_15h[,-252] <- actD_15h_t

actD_0hr$means <- apply(actD_0hr[-252],1,function(x){weighted.mean(x=8:250,w=x[-c(1:8)])})
actD_1hr$means <- apply(actD_1hr[-252],1,function(x){weighted.mean(x=8:250,w=x[-c(1:8)])})
actD_3hr$means <- apply(actD_3hr[-252],1,function(x){weighted.mean(x=8:250,w=x[-c(1:8)])})
actD_7hr$means <- apply(actD_7hr[-252],1,function(x){weighted.mean(x=8:250,w=x[-c(1:8)])})
actD_15h$means <- apply(actD_15h[-252],1,function(x){weighted.mean(x=8:250,w=x[-c(1:8)])})


all <- join_all(list(actD_0hr,
	actD_1hr,
	actD_3hr,
	actD_7hr,
	actD_15h),by="accession")
all <- all[complete.cases(all),]
print(nrow(all))

actD_0hr$tp = "tp0"
actD_1hr$tp = "tp1"
actD_3hr$tp = "tp3"
actD_7hr$tp = "tp7"
actD_15h$tp = "t15"

actD_0hr <- actD_0hr[which(actD_0hr$accession %in% all$accession),]
actD_1hr <- actD_1hr[which(actD_1hr$accession %in% all$accession),]
actD_3hr <- actD_3hr[which(actD_3hr$accession %in% all$accession),]
actD_7hr <- actD_7hr[which(actD_7hr$accession %in% all$accession),]
actD_15h <- actD_15h[which(actD_15h$accession %in% all$accession),]

print(median(actD_0hr$means))
print(median(actD_1hr$means))
print(median(actD_3hr$means))
print(median(actD_7hr$means))
print(median(actD_15h$means))

mn <- rbind(
	actD_0hr[,c(253,254)],
	actD_1hr[,c(253,254)],
	actD_3hr[,c(253,254)],
	actD_7hr[,c(253,254)],
	actD_15h[,c(253,254)])

m$tail_length <- c(4:124,125)*2
all_dat_melt <- gather(m,tp,val,-tail_length)
# all_dat_melt$norm <- rep(c(0.000585759,0.000385634,0.000497605,0.001027441,0.002286478),each=nrow(all_dat_melt)/5)
all_dat_melt$norm <- rep(1/c(1243.8548, 1184.7926, 1058.5599,  401.1526,  208.4624),each=nrow(all_dat_melt)/5) #from the RNAseq data, as in 5E.
all_dat_melt_tb <- as_tibble(all_dat_melt)
all_dat_melt_tb$normabr <- all_dat_melt_tb$val/all_dat_melt_tb$norm/sum(all_dat_melt_tb[all_dat_melt_tb$tp == "tp0",]$val/all_dat_melt_tb[all_dat_melt_tb$tp == "tp0",]$norm)
# all_dat_melt_tb <- group_by(all_dat_melt_tb, tp) %>% mutate(normabr = normab/sum(normab,na.rm = TRUE)) %>% ungroup()

p1<-ggplot(all_dat_melt_tb[!all_dat_melt$tail_length==250,],aes(x=tail_length,y=normabr,color=tp)) + 
geom_step(size=0.5) + 
scale_y_continuous(name=NULL,expand=c(0,0),limits=c(0,NA)) + 
theme_tim() + 
scale_x_continuous(name=NULL,limits=c(8,250), labels = c(8, 1:5 * 50), breaks = c(8, 1:5 * 50)) + 
scale_colour_manual(values=cols)+
theme(axis.line.x = element_blank()) + 
geom_segment(aes(x = 8, xend = 250, y = 0, yend = 0),color = "black")

print(all_dat_melt_tb[all_dat_melt_tb$tail_length == 250,])
p2<-ggplot(all_dat_melt[!all_dat_melt$tail_length==250,],aes(x=tail_length,y=val,color=tp)) + 
geom_step(size=0.5) + 
scale_y_continuous(name=NULL,expand=c(0,0),limits=c(0,NA)) + 
theme_tim() + 
scale_x_continuous(name=NULL,limits=c(8,250), labels = c(8, 1:5 * 50), breaks = c(8, 1:5 * 50)) + 
scale_colour_manual(values=cols) +
theme(axis.line.x = element_blank()) + 
geom_segment(aes(x = 8, xend = 250, y = 0, yend = 0),color = "black")

print(head(mn))

p3<-ggplot(mn,aes(x=means,y=..density..,colour=tp)) + 
stat_bin(geom="step",bins=118,position='identity',size=0.5) + 
theme_tim() + 
scale_x_continuous(name=NULL,limits=c(8,250), labels = c(8, 1:5 * 50), breaks = c(8, 1:5 * 50)) + 
scale_y_continuous(name=NULL,expand=c(0,0),limits=c(0,NA)) + 
scale_colour_manual(values=cols)+
theme(axis.line.x = element_blank()) + 
geom_segment(aes(x = 8, xend = 250, y = 0, yend = 0),color = "black")

ggsave(plot=p1,file="figures/version_9/Fig5Bpanel1.pdf",width=2,height=1)
ggsave(plot=p2,file="figures/version_9/Fig5Bpanel2.pdf",width=2,height=1)
ggsave(plot=p3,file="figures/version_9/Fig5C.pdf",width=2,height=1)
