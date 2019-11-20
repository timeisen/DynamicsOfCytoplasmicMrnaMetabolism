#Read data in 
#SWE uses 
library(ggplot2)
library(tictoc)
args <- commandArgs(trailingOnly=TRUE)
# sc <- read.table("/lab/solexa_bartel/teisen/RNAseq/kinetics_of_translation/final_analyses/analysis/processed_files/miR-155_RNA_minus_tags_10rpm_cutoff_with_stds_scaled.txt",head=TRUE,sep="\t")
options(warn=1)
sc <- read.table(args[1],head=TRUE,sep="\t")
# sc <- sc[,c(1,8:13)]
colnames(sc)<-c("accession","x40m","x1h","x2h","x4h","xSS")
stds<-read.table("/lab/solexa_bartel/teisen/RNAseq/Annotation_files/TAIL_seq/list_of_standards_SWE_notation.txt")

sc <- sc[which(!sc$accession %in% stds$V1),]

hl<-read.table("/lab/solexa_bartel/eichhorn/5EU_RPF/miR-1_miR-155_analysis/Halflife_analysis/Halflife_measurements/miR-155_RNA_minusdox_halflife_SS.txt",head=TRUE,sep="\t")

# sc <- head(sc,50)
objective<-function(params,time,data,offset){
	alpha = exp(params[1])
	beta = exp(params[2])
    model = sapply(time,function(t,alpha,beta,offset){alpha/beta*(1-exp(-beta*(t-offset)))},alpha=alpha,beta=beta,offset=offset)
	sqs<-sum((log(model)-log(data))^2)
	return(sqs)
}

optimize<-function(time,data,guess,offset){
	test_param <- log(c(guess,1))
	solve<-optim(
		par=test_param,
		fn=objective,
		method="L-BFGS-B",
		lower = log(rep(1E-15,2)),
		upper = log(c(rep(1E15,2))),
		data=data[-c(5)],
		time=time,
		offset = offset)
		return(c(exp(solve$par),solve$value))}

time_pts<-c(40,60,120,240) #it seems like a good offset is 56minutes
# fits<-data.frame(t(apply(sc[,-1],1,function(x){optimize(time_pts,x)})))
guess = median(sc[,2])
allfits <- NULL
rs <- NULL
for (offset in 1:39){
	fits <- data.frame(t(apply(sc[,-1],1,function(x,time_pts,guess,offset){optimize(time_pts,x,guess,offset)},time=time_pts,guess=guess,offset=offset)))
	fits <- cbind(sc[1],fits)
	colnames(fits)<-c("accession","alpha_t","beta_t","residual_t")
	fits$halflife_t <- log(2)/fits$beta_t/60
	hlo <- merge(hl,fits,by="accession")
	fits$offset <- offset
	allfits <- rbind(allfits,fits)
	rs <- rbind(rs,c(sum(fits$residual_t),offset))
}
toff <- rs[which.min(rs[,1]),2]
fits <- allfits[which(allfits$offset == toff),]
print(fits)
print(rs[toff,1])
write.table(fits,file=args[2],quote=FALSE,sep="\t",row.names=FALSE)

# p1 <- ggplot(hl,aes(x=halflife,y=halflife_t,color=log(residual_t)))+scale_x_continuous(trans="log10",limits=c(0.001,100))+scale_y_continuous(trans="log10",limits=c(0.001,100))+geom_point()+coord_fixed()
# print(cor(hl$halflife,hl$halflife_t,method="spearman"))
# print(cor(hl$alpha,hl$alpha_t,method="spearman"))

#Scale the steady state value7.58771929657