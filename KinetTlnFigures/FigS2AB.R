library(ggplot2)
source("/lab/solexa_bartel/teisen/RNAseq/Scripts/general/ggplot_theme.R")


stds <- read.table("/lab/solexa_bartel/teisen/RNAseq/Annotation_files/TAIL_seq/list_of_standards_SWE_notation.txt")
std_set_1 <- stds[c(1,3,5,6,8,9,11),]
std_set_2 <- stds[c(2,4,7,10),]

lis_data<-NULL

#Stephen's data
lis_data[[1]] <- read.table("/lab/solexa_bartel/eichhorn/5EU_RPF/miR-1_miR-155_timecourse_tail_analysis/unfiltered_analysis/Tail_lengths/TAGTGC-1_mediantails.txt",head=TRUE,sep="\t")
lis_data[[2]] <- read.table("/lab/solexa_bartel/eichhorn/5EU_RPF/miR-1_miR-155_timecourse_tail_analysis/unfiltered_analysis/Tail_lengths/GCTACA-1_mediantails.txt",head=TRUE,sep="\t")
lis_data[[3]] <- read.table("/lab/solexa_bartel/eichhorn/5EU_RPF/miR-1_miR-155_timecourse_tail_analysis/unfiltered_analysis/Tail_lengths/AATCCG-1_mediantails.txt",head=TRUE,sep="\t")
lis_data[[4]] <- read.table("/lab/solexa_bartel/eichhorn/5EU_RPF/miR-1_miR-155_timecourse_tail_analysis/unfiltered_analysis/Tail_lengths/ACTGGT-1_mediantails.txt",head=TRUE,sep="\t")
lis_data[[5]] <- read.table("/lab/solexa_bartel/eichhorn/5EU_RPF/miR-1_miR-155_timecourse_tail_analysis/unfiltered_analysis/Tail_lengths/CGGTTA-1_mediantails.txt",head=TRUE,sep="\t")
lis_data[[6]] <- read.table("/lab/solexa_bartel/eichhorn/5EU_RPF/miR-1_miR-155_timecourse_tail_analysis/unfiltered_analysis/Tail_lengths/GTCTAG-1_mediantails.txt",head=TRUE,sep="\t")
lis_data[[7]] <- read.table("/lab/solexa_bartel/eichhorn/5EU_RPF/miR-1_miR-155_timecourse_tail_analysis/unfiltered_analysis/Tail_lengths/TGAACT-1_mediantails.txt",head=TRUE,sep="\t")
lis_data[[8]] <- read.table("/lab/solexa_bartel/eichhorn/5EU_RPF/miR-1_miR-155_timecourse_tail_analysis/unfiltered_analysis/Tail_lengths/CTAGTC-1_mediantails.txt",head=TRUE,sep="\t")
lis_data[[9]] <- read.table("/lab/solexa_bartel/eichhorn/5EU_RPF/miR-1_miR-155_timecourse_tail_analysis/unfiltered_analysis/Tail_lengths/ATCGAA-1_mediantails.txt",head=TRUE,sep="\t")
lis_data[[10]] <- read.table("/lab/solexa_bartel/eichhorn/5EU_RPF/miR-1_miR-155_timecourse_tail_analysis/unfiltered_analysis/Tail_lengths/CAGCGT-1_mediantails.txt",head=TRUE,sep="\t")

#Tim's data
lis_data[[11]] <- read.table("/lab/solexa_bartel/eichhorn/5EU_RPF/miR-1_miR-155_timecourse_tail_analysis/unfiltered_analysis/Tail_lengths/TAGTGC-2_mediantails.txt",head=TRUE,sep="\t")
lis_data[[12]] <- read.table("/lab/solexa_bartel/eichhorn/5EU_RPF/miR-1_miR-155_timecourse_tail_analysis/unfiltered_analysis/Tail_lengths/GCTACA-2_mediantails.txt",head=TRUE,sep="\t")
lis_data[[13]] <- read.table("/lab/solexa_bartel/eichhorn/5EU_RPF/miR-1_miR-155_timecourse_tail_analysis/unfiltered_analysis/Tail_lengths/AATCCG-2_mediantails.txt",head=TRUE,sep="\t")
lis_data[[14]] <- read.table("/lab/solexa_bartel/eichhorn/5EU_RPF/miR-1_miR-155_timecourse_tail_analysis/unfiltered_analysis/Tail_lengths/ACTGGT-2_mediantails.txt",head=TRUE,sep="\t")
lis_data[[15]] <- read.table("/lab/solexa_bartel/eichhorn/5EU_RPF/miR-1_miR-155_timecourse_tail_analysis/unfiltered_analysis/Tail_lengths/TGTCAC-2_mediantails.txt",head=TRUE,sep="\t")
lis_data[[16]] <- read.table("/lab/solexa_bartel/eichhorn/5EU_RPF/miR-1_miR-155_timecourse_tail_analysis/unfiltered_analysis/Tail_lengths/CGGTTA-2_mediantails.txt",head=TRUE,sep="\t")
lis_data[[17]] <- read.table("/lab/solexa_bartel/eichhorn/5EU_RPF/miR-1_miR-155_timecourse_tail_analysis/unfiltered_analysis/Tail_lengths/GTCTAG-2_mediantails.txt",head=TRUE,sep="\t")
lis_data[[18]] <- read.table("/lab/solexa_bartel/eichhorn/5EU_RPF/miR-1_miR-155_timecourse_tail_analysis/unfiltered_analysis/Tail_lengths/TGAACT-2_mediantails.txt",head=TRUE,sep="\t")
lis_data[[19]] <- read.table("/lab/solexa_bartel/eichhorn/5EU_RPF/miR-1_miR-155_timecourse_tail_analysis/unfiltered_analysis/Tail_lengths/CTAGTC-2_mediantails.txt",head=TRUE,sep="\t")
lis_data[[20]] <- read.table("/lab/solexa_bartel/eichhorn/5EU_RPF/miR-1_miR-155_timecourse_tail_analysis/unfiltered_analysis/Tail_lengths/ATCGAA-2_mediantails.txt",head=TRUE,sep="\t")
lis_data[[21]] <- read.table("/lab/solexa_bartel/eichhorn/5EU_RPF/miR-1_miR-155_timecourse_tail_analysis/unfiltered_analysis/Tail_lengths/GCAATT-2_mediantails.txt",head=TRUE,sep="\t")
lis_data[[22]] <- read.table("/lab/solexa_bartel/eichhorn/5EU_RPF/miR-1_miR-155_timecourse_tail_analysis/unfiltered_analysis/Tail_lengths/CAGCGT-2_mediantails.txt",head=TRUE,sep="\t")

#ActD data
lis_data[[23]] <- read.table(Sys.glob(file.path("/lab/solexa_bartel/teisen/Tail-seq/miR-1_actD_samples/data_files_20171030/0hr-dox_input_17_10_31/*mediantails.txt")),head=TRUE,sep="\t")
lis_data[[24]] <- read.table(Sys.glob(file.path("/lab/solexa_bartel/teisen/Tail-seq/miR-1_actD_samples/data_files_20171030/1hr-dox_input_17_11_03/*mediantails.txt")),head=TRUE,sep="\t")
lis_data[[25]] <- read.table(Sys.glob(file.path("/lab/solexa_bartel/teisen/Tail-seq/miR-1_actD_samples/data_files_20171030/3hr-dox_input_17_10_30/*mediantails.txt")),head=TRUE,sep="\t")
lis_data[[26]] <- read.table(Sys.glob(file.path("/lab/solexa_bartel/teisen/Tail-seq/miR-1_actD_samples/data_files_20171030/7hr-dox_input_17_11_03/*mediantails.txt")),head=TRUE,sep="\t")
lis_data[[27]] <- read.table(Sys.glob(file.path("/lab/solexa_bartel/teisen/Tail-seq/miR-1_actD_samples/data_files_20171030/15hr-dox_input_17_11_03/*mediantails.txt")),head=TRUE,sep="\t")
# lis_data[[28]] <- read.table(Sys.glob(file.path("/lab/solexa_bartel/teisen/Tail-seq/miR-1_actD_samples/data_files_20171030/0hr+dox_input_17_10_30/*mediantails.txt")),head=TRUE,sep="\t")
# lis_data[[29]] <- read.table(Sys.glob(file.path("/lab/solexa_bartel/teisen/Tail-seq/miR-1_actD_samples/data_files_20171030/1hr+dox_input_17_10_30/*mediantails.txt")),head=TRUE,sep="\t")
# lis_data[[30]] <- read.table(Sys.glob(file.path("/lab/solexa_bartel/teisen/Tail-seq/miR-1_actD_samples/data_files_20171030/3hr+dox_input_17_10_30/*mediantails.txt")),head=TRUE,sep="\t")
# lis_data[[31]] <- read.table(Sys.glob(file.path("/lab/solexa_bartel/teisen/Tail-seq/miR-1_actD_samples/data_files_20171030/7hr+dox_input_17_10_31/*mediantails.txt")),head=TRUE,sep="\t")
# lis_data[[32]] <- read.table(Sys.glob(file.path("/lab/solexa_bartel/teisen/Tail-seq/miR-1_actD_samples/data_files_20171030/15hr+dox_input_17_10_31/*mediantails.txt")),head=TRUE,sep="\t")
lis_data[[28]] <- read.table(Sys.glob(file.path("/lab/solexa_bartel/teisen/Tail-seq/miR-1_actD_samples/data_files_20171030/0hr-dox_eluate_17_10_30/*mediantails.txt")),head=TRUE,sep="\t")
lis_data[[29]] <- read.table(Sys.glob(file.path("/lab/solexa_bartel/teisen/Tail-seq/miR-1_actD_samples/data_files_20171030/1hr-dox_eluate_17_10_30/*mediantails.txt")),head=TRUE,sep="\t")
lis_data[[30]] <- read.table(Sys.glob(file.path("/lab/solexa_bartel/teisen/Tail-seq/miR-1_actD_samples/data_files_20171030/3hr-dox_eluate_17_10_30/*mediantails.txt")),head=TRUE,sep="\t")
lis_data[[31]] <- read.table(Sys.glob(file.path("/lab/solexa_bartel/teisen/Tail-seq/miR-1_actD_samples/data_files_20171030/7hr-dox_eluate_17_10_30/*mediantails.txt")),head=TRUE,sep="\t")
lis_data[[32]] <- read.table(Sys.glob(file.path("/lab/solexa_bartel/teisen/Tail-seq/miR-1_actD_samples/data_files_20171030/15hr-dox_eluate_17_10_30/*mediantails.txt")),head=TRUE,sep="\t")
lis_data[[33]] <- read.table(Sys.glob(file.path("/lab/solexa_bartel/teisen/Tail-seq/miR-1_actD_samples/data_files_20171030/total-dox_eluate_17_10_30/*mediantails.txt")),head=TRUE,sep="\t")
# lis_data[[39]] <- read.table(Sys.glob(file.path("/lab/solexa_bartel/teisen/Tail-seq/miR-1_actD_samples/data_files_20171030/0hr+dox_eluate_17_10_30/*mediantails.txt")),head=TRUE,sep="\t")
# lis_data[[40]] <- read.table(Sys.glob(file.path("/lab/solexa_bartel/teisen/Tail-seq/miR-1_actD_samples/data_files_20171030/1hr+dox_eluate_17_10_30/*mediantails.txt")),head=TRUE,sep="\t")
# lis_data[[41]] <- read.table(Sys.glob(file.path("/lab/solexa_bartel/teisen/Tail-seq/miR-1_actD_samples/data_files_20171030/3hr+dox_eluate_17_10_31/*mediantails.txt")),head=TRUE,sep="\t")
# lis_data[[42]] <- read.table(Sys.glob(file.path("/lab/solexa_bartel/teisen/Tail-seq/miR-1_actD_samples/data_files_20171030/7hr+dox_eluate_17_10_31/*mediantails.txt")),head=TRUE,sep="\t")
# lis_data[[43]] <- read.table(Sys.glob(file.path("/lab/solexa_bartel/teisen/Tail-seq/miR-1_actD_samples/data_files_20171030/15hr+dox_eluate_17_10_31/*mediantails.txt")),head=TRUE,sep="\t")
# lis_data[[44]] <- read.table(Sys.glob(file.path("/lab/solexa_bartel/teisen/Tail-seq/miR-1_actD_samples/data_files_20171030/total+dox_eluate_17_10_30/*mediantails.txt")),head=TRUE,sep="\t")

#Nuclear data
# lis_data[[45]] <- read.table("/lab/solexa_bartel/teisen/Tail-seq/nuclear_tails/data/total_17_12_21/tail_lengthsTAGTGC-2_mediantails.txt",head=TRUE,sep="\t")
# lis_data[[46]] <- read.table("/lab/solexa_bartel/teisen/Tail-seq/nuclear_tails/data/chromatin_17_12_21/tail_lengthsTGTCAC-2_mediantails.txt",head=TRUE,sep="\t")
# lis_data[[47]] <- read.table("/lab/solexa_bartel/teisen/Tail-seq/nuclear_tails/data/cytoplasm_17_12_21/tail_lengthsGCTACA-2_mediantails.txt",head=TRUE,sep="\t")

arr_std_set_1 <- NULL
arr_std_set_2 <- NULL

arr_std_expr <- data.frame(
	Accession.number = c("NM_A10_7",
					 "NM_A50_7",
					 "NM_A100_7",
					 "NM_A150_7",
					 "NM_A200_7",
					 "NM_A250_7",
					 "NM_A300_7",
					 "NM_A10_6",
					 "NM_A30_6",
					 "NM_A110_6",
					 "NM_A210_6"),
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
						1.135))

for(i in 1:length(lis_data)){
	
	temp_set_1 <- lis_data[[i]][which(lis_data[[i]]$Accession.number %in% std_set_1),]
	temp_set_2 <- lis_data[[i]][which(lis_data[[i]]$Accession.number %in% std_set_2),]

	#normalize the data to the A10 and to the scaled_abundance
	temp_set_1$Tags <- temp_set_1$Tags/temp_set_1[which(temp_set_1$Accession.number == "NM_A10_7"),]$Tags
	temp_set_2$Tags <- temp_set_2$Tags/temp_set_2[which(temp_set_2$Accession.number == "NM_A10_6"),]$Tags

	temp_set_1 <- merge(temp_set_1,arr_std_expr,by="Accession.number")
	temp_set_2 <- merge(temp_set_2,arr_std_expr,by="Accession.number")

	temp_set_1$Tags <- temp_set_1$Tags/temp_set_1$scaled_abundance
	temp_set_2$Tags <- temp_set_2$Tags/temp_set_2$scaled_abundance

	temp_set_1$scaled_abundance <- NULL
	temp_set_2$scaled_abundance <- NULL

	#store data
	arr_std_set_1 <- rbind(
		temp_set_1,
		arr_std_set_1)
	arr_std_set_2 <- rbind(
		temp_set_2,
		arr_std_set_2)

	}


arr_std_set_1$set <- "Set 1"
arr_std_set_2$set <- "Set 2"

arr_std_set_1$Accession.number <- factor(arr_std_set_1$Accession.number,levels=c(
																	"NM_A10_7",
																	"NM_A50_7",
																	"NM_A100_7",
																	"NM_A150_7",
																	"NM_A200_7",
																	"NM_A250_7",
																	"NM_A300_7"))

arr_std_set_2$Accession.number <- factor(arr_std_set_2$Accession.number,levels=c(
																	"NM_A10_6",
																	"NM_A30_6",
																	"NM_A110_6",
																	"NM_A210_6"))

arr_std_values <- data.frame(
	Accession.number = c("NM_A10_7",
						 "NM_A50_7",
						 "NM_A100_7",
						 "NM_A150_7",
						 "NM_A200_7",
						 "NM_A250_7",
						 "NM_A300_7",
						 "NM_A10_6",
						 "NM_A30_6",
						 "NM_A110_6",
						 "NM_A210_6"),
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
					210),
	set = c(rep("Set 1",7),rep("Set 2",4))
	)

arr_stds <- rbind(arr_std_set_1,arr_std_set_2)
print(arr_stds)
p1<-ggplot(arr_stds,aes(x=Accession.number,y=Mean.length),shape=1,color="black")+geom_point(size=0.5)+theme_tim()+
	facet_grid(.~set,scales="free_x",space="free_x")+geom_point(data = arr_std_values,aes(x=Accession.number,y=Mean.length),shape=3,color="red",size = 2)+
	scale_y_continuous(limits=c(0,350),name=NULL)+scale_x_discrete(labels=c("10","57","107","160","215","273","324"))

p2<-ggplot(arr_stds[which(arr_stds$Accession.number != "NM_A30_6"),],aes(x=Accession.number,y=Tags),shape=1,color="black",size=0.5)+
	geom_point(size=0.5,position=position_dodge(width = 0.3))+
	theme_tim()+
	facet_grid(.~set,scales="free_x",space="free_x")+
	scale_y_continuous(limits=c(0.01,10),labels=log_ticks(0.01,10)[[1]],breaks=log_ticks(0.01,10)[[2]],name=NULL,trans="log10")+scale_x_discrete(labels=c("10","57","107","160","215","273","324"))

ggsave(plot=p2,file="/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/figures/version_6/FigS2A.pdf",useDingbats=FALSE,width=3,height=2)
ggsave(plot=p1,file="/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/figures/version_6/FigS2B.pdf",useDingbats=FALSE,width=3,height=2)
