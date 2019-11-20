#reads in raw count data, generates new datasets by sampling with replacement.
args <- commandArgs(trailingOnly=TRUE)
dat <- read.table(args[1],sep="\t",head = TRUE)
print("does the input dataset have a header?")
outdir <- args[2]
bn <- tools::file_path_sans_ext(basename(args[1]))

bootstrap <- function(df_m){
	vec_m <- c(as.matrix(df_m[,-1]))
	int_s = sum(vec_m)
	mat_m_sampled <- matrix(rmultinom(n = 1,size =  int_s, prob = vec_m),nrow=nrow(df_m),ncol=ncol(df_m)-1)
	df_m_sampled <- cbind(df_m[,1],data.frame(mat_m_sampled))
	return(df_m_sampled)}

for(i in 1:10){
	df_m_sampled <- bootstrap(dat)
	filename = paste0(outdir,bn,"_BOOTSTRAP_",i,".txt")
	write.table(df_m_sampled,file=filename,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)}



#tp40m<-read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/single_tag_files/40min_minus_tag_reformat_no_cutoff_miR-155.txt",sep="\t")

#/lab/solexa_bartel/teisen/Tail-seq/PalSeqDirectLig/Data/Run2AllDataReProcess/CombinedData/miR-155_uninduced/tail_lengths_annot_end_only_reformat.txt