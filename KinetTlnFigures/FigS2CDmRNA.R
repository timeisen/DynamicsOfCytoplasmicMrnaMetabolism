library(ggplot2)
library(plyr)
library(reshape2)
library(cowplot)
options("scipen"=5)
source("/lab/solexa_bartel/teisen/RNAseq/Scripts/general/ggplot_theme.R")

PAL_V1_1__ <-read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/other_analyses/directLig/process_Stds/STDS_PAL-seq_miR-1_minus_2hr_compiled.txt",nrow=14)
PAL_V1_155 <-read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/other_analyses/directLig/process_Stds/STDS_PAL-seq_miR-155_minus_2hr_compiled.txt",nrow=14)
TAIL___1__ <-read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/other_analyses/directLig/process_Stds/STDS_2hr_minus_v8_pcr_miR1_v2_compiled.txt")
TAIL___155 <-read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/other_analyses/directLig/process_Stds/STDS_2hr_minus_v8_pcr_v2_compiled.txt")
PAL_SS_1__ <-read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/other_analyses/directLig/process_Stds/miR1_uninduced_R0_compiled20190715.txt") 
PAL_SS_155 <-read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/other_analyses/directLig/process_Stds/miR155_uninduced_R0_compiled20190715.txt") 


colnames(PAL_V1_1__) <- c("accession","mean","abundance")
colnames(PAL_V1_155) <- c("accession","mean","abundance")
colnames(TAIL___1__) <- c("accession","mean","abundance")
colnames(TAIL___155) <- c("accession","mean","abundance")
colnames(PAL_SS_1__) <- c("accession","mean","abundance")
colnames(PAL_SS_155) <- c("accession","mean","abundance")


palAcc <- c(
	"STD10",
	"STD100",
	"A100_7",
	"A10_7",
	"STD150",
	"A150_7",
	"STD200",
	"A200_7",
	"STD250",
	"A250_7",
	"STD300",
	"A300_7",
	"STD50",
	"A50_7")

PAL_V1_1__$accession = palAcc
PAL_V1_155$accession = palAcc


PAL_V1_1__$abundance = PAL_V1_1__$abundance/PAL_V1_1__[PAL_V1_1__$accession=="A100_7",]$abundance
PAL_V1_155$abundance = PAL_V1_155$abundance/PAL_V1_155[PAL_V1_155$accession=="A100_7",]$abundance
TAIL___1__$abundance = TAIL___1__$abundance/TAIL___1__[TAIL___1__$accession=="A100_7",]$abundance
TAIL___155$abundance = TAIL___155$abundance/TAIL___155[TAIL___155$accession=="A100_7",]$abundance
PAL_SS_1__$abundance = PAL_SS_1__$abundance/PAL_SS_1__[PAL_SS_1__$accession=="A100_7",]$abundance
PAL_SS_155$abundance = PAL_SS_155$abundance/PAL_SS_155[PAL_SS_155$accession=="A100_7",]$abundance


PAL_V1_1__$Data = "PAL_V1_1__"
TAIL___1__$Data = "TAIL___1__"
PAL_V1_155$Data = "PAL_V1_155"
TAIL___155$Data = "TAIL___155"
PAL_SS_1__$Data = "PAL_SS_1__"
PAL_SS_155$Data = "PAL_SS_155"



arr_stds <- rbind(
	PAL_V1_1__,
	TAIL___1__,
	PAL_V1_155,
	TAIL___155,
	PAL_SS_1__,
	PAL_SS_155)

print(dim(arr_stds))


arr_std_expr <- data.frame(
	accession = c(   "A10_7",
					 "A50_7",
					 "A100_7",
					 "A150_7",
					 "A200_7",
					 "A250_7",
					 "A300_7",
					 "A10_6",
					 "A30_6",
					 "A110_6",
					 "A210_6",
					 "STD10",
					 "STD50",
					 "STD100",
					 "STD150",
					 "STD200",
					 "STD250",
					 "STD300"
					 ),
	scaled_abundance = c(1,
						1.268,
						1.328,
						1.346,
						1.040,
						0.147,
						0.0951,
						1,
						2.729,
						0.259,
						1.135,
						1,
						1.416,
						1.933,
						0.943,
						0.897,
						0.249,
						0.203)/1.328)
print(dim(arr_std_expr))
arr_stds <- merge(arr_stds,arr_std_expr,by="accession")
print(dim(arr_stds))

arr_stds$Tags = arr_stds$abundance/arr_stds$scaled_abundance

arr_std_values <- data.frame(
	accession = c(
						 "A10_7",
						 "A50_7",
						 "A100_7",
						 "A150_7",
						 "A200_7",
						 "A250_7",
						 "A300_7",
						 "A10_6",
						 "A30_6",
						 "A110_6",
						 "A210_6",
						 "STD10",
					     "STD50",
					     "STD100",
					     "STD150",
					     "STD200",
					     "STD250",
					     "STD300"),
	Mean.length = c(10,
					57,
					107,
					160,
					215,
					273,
					324,
					10,
					30,
					110,
					210,
					10,
					56,
					107,
					167,
					205,
					263,
					302),
	set = c(rep("Set 1",7),rep("Set 2",4),rep("Set 3",7))
	)
arr_stds <- merge(arr_stds,arr_std_values,by="accession")
print(dim(arr_stds))

set1 <- c(
																	"A10_7",
																	"A50_7",
																	"A100_7",
																	"A150_7",
																	"A200_7",
																	"A250_7",
																	"A300_7")

set2 <- c(
																	"A10_6",
																	"A30_6",
																	"A110_6",
																	"A210_6")

set3 <- c(
																	"STD10",
																	"STD50",
																	"STD100",
																	"STD150",
																	"STD200",
																	"STD250",
																	"STD300")




arr_std_set_1<- arr_stds[which(arr_stds$accession %in% set1),]
arr_std_set_2<- arr_stds[which(arr_stds$accession %in% set2),]
arr_std_set_3<- arr_stds[which(arr_stds$accession %in% set3),]

arr_std_set_1$set <- "Set 1"
arr_std_set_2$set <- "Set 2"
arr_std_set_3$set <- "Set 3"


arr_std_set_1$accession <- factor(arr_std_set_1$accession,levels=c(
																	"A10_7",
																	"A50_7",
																	"A100_7",
																	"A150_7",
																	"A200_7",
																	"A250_7",
																	"A300_7"))

arr_std_set_2$accession <- factor(arr_std_set_2$accession,levels=c(
																	"A10_6",
																	"A30_6",
																	"A110_6",
																	"A210_6"))

arr_std_set_3$accession <- factor(arr_std_set_3$accession,levels=c(
																	"STD10",
																	"STD50",
																	"STD100",
																	"STD150",
																	"STD200",
																	"STD250",
																	"STD300"))


arr_stds <- rbind(arr_std_set_1,arr_std_set_2,arr_std_set_3)

arr_std_values$accession <- factor(arr_std_values$accession,levels=c(
																	"A10_7",
																	"A50_7",
																	"A100_7",
																	"A150_7",
																	"A200_7",
																	"A250_7",
																	"A300_7",
																	"A10_6",
																	"A30_6",
																	"A110_6",
																	"A210_6",
																	"STD10",
																	"STD50",
																	"STD100",
																	"STD150",
																	"STD200",
																	"STD250",
																	"STD300"))
print(arr_stds)
p1<-ggplot(arr_stds,aes(x=accession,y=mean,color=Data),shape=1)+
	geom_point(size=0.5,position=position_dodge(width = 0.3))+
	theme_tim()+
	facet_grid(.~set,scales="free_x",space="free_x")+
	geom_point(data = arr_std_values,aes(x=accession,y=Mean.length),shape=3,color="red",size = 2)+
	scale_y_continuous(limits=c(0,350),name=NULL)+scale_x_discrete("abundance",labels=c("10","57","107","160","215","273","324"))

print("Set X axis labels manually for these plots.")

p2<-ggplot(arr_stds[which(arr_stds$accession != "A30_6"),],aes(x=accession,y=Tags,color=Data),shape=1,color="black",size=0.5)+
	geom_point(size=0.5,position=position_dodge(width = 0.3))+
	theme_tim()+
	facet_grid(.~set,scales="free_x",space="free_x")+
	scale_y_continuous(limits=c(0.01,100),labels=log_ticks(0.01,100)[[1]],breaks=log_ticks(0.01,100)[[2]],name=NULL,trans="log10",expand=c(0,0))+scale_x_discrete(labels=c("10","57","107","160","215","273","324"))

ggsave(plot=p1,file="/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/figures/version_9/FigS5E_v2.pdf",useDingbats=FALSE,width=3.25,height=2)
ggsave(plot=p2,file="/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/figures/version_9/FigS5D_v2.pdf",useDingbats=FALSE,width=3.25,height=2)
ggsave(plot=get_legend(p1 + theme(legend.position = "right")),file="/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/figures/version_9/FigS5E_v2_legend.pdf",useDingbats=FALSE,width=3.25,height=2)
ggsave(plot=get_legend(p2 + theme(legend.position = "right")),file="/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/figures/version_9/FigS5D_v2_legend.pdf",useDingbats=FALSE,width=3.25,height=2)
