#Read data in 
library(ggplot2)
library(tictoc)
args <- commandArgs(trailingOnly=TRUE)
# sc <- read.table("/lab/solexa_bartel/teisen/RNAseq/kinetics_of_translation/final_analyses/analysis/processed_files/miR-155_RNA_minus_tags_10rpm_cutoff_with_stds_scaled.txt",head=TRUE,sep="\t")
sc <- read.table(args[1],head=TRUE,sep="\t")
options(warn=1)

# sc <- sc[,c(1,8:13)]
colnames(sc)<-c("accession","x40m","x1h","x2h","x4h","xSS")
stds<-read.table("/lab/solexa_bartel/teisen/RNAseq/Annotation_files/TAIL_seq/list_of_standards_SWE_notation.txt")

sc <- sc[which(!sc$accession %in% stds$V1),]

hl<-read.table("/lab/solexa_bartel/eichhorn/5EU_RPF/miR-1_miR-155_analysis/Halflife_analysis/Halflife_measurements/miR-155_RNA_minusdox_halflife_SS.txt",head=TRUE,sep="\t")

# sc <- head(sc,50)
objective<-function(params,time,data,gene_no){
	alpha = exp(params[1:gene_no])
	beta = exp(params[(gene_no+1):(gene_no*2)])
	offset = params[gene_no*2+1]
    model = sapply(time,function(t,alpha,beta,offset){alpha/beta*(1-exp(-beta*(t-offset)))},alpha=alpha,beta=beta,offset=offset)
	sqs<-sum((log(model)-log(data))^2)
	return(sqs)
}
###Analytic gradient for optimization efficiency
grAnalytic <- function(params,time,data,gene_no){
	alpha = exp(params[1:gene_no]) #because these are exponentiated here, I need to multiply the final "alphas" by alpha, same with beta.
	beta = exp(params[(gene_no+1):(gene_no*2)])
	offset = params[gene_no*2+1]
    model = sapply(time,function(t,alpha,beta,offset){alpha/beta*(1-exp(-beta*(t-offset)))},alpha=alpha,beta=beta,offset=offset)
	preFactor = 2*(log(model)-log(data))/model
	#derivative from alphas
	alphas = rowSums(preFactor*sapply(time,function(t,alpha,beta,offset){(1-exp(-beta*(t-offset)))/beta},alpha=alpha,beta=beta,offset=offset))*alpha
	#derivative of betas
	betas = rowSums(preFactor*sapply(time,function(t,alpha,beta,offset){-alpha*(1-exp(-beta*(t-offset)))/beta^2-alpha*exp(-beta*(t-offset))*(offset-t)/beta},alpha=alpha,beta=beta,offset=offset))*beta
	#derivative of offset
	offset = sum(preFactor*sapply(time,function(t,alpha,beta,offset){-alpha*exp(-beta*(t-offset))},alpha=alpha,beta=beta,offset=offset))
	gradient = c(alphas,betas,offset)
    return(gradient)
}

optimize<-function(time,data,guess){
	gene_no = nrow(data)
	test_param <- c(log(rep(c(guess,1),each = gene_no)),32)
	for(i in 1:2){
	solve<-optim(
		par=test_param,
		fn=objective,
		method="L-BFGS-B",
		gr=grAnalytic,
		lower = c(log(rep(c(1E-10,1E-5),each = gene_no)),1),
		upper = c(log(rep(c(1E2,1E2),each = gene_no)),39),
		data=data[-c(1)],
		time=time,
		gene_no = gene_no)
	test_param = solve$par}
	return(c(exp(solve$par[-length(solve$par)]),solve$par[length(solve$par)],solve$value))}

time_pts<-c(40,60,120,240,6000) #it seems like a good offset is 56minutes
guess = median(sc[,2])
pars <- optimize(time_pts,sc,guess)
fits <- data.frame(accession = sc[,1],alpha_t = pars[1:nrow(sc)],beta_t = pars[(nrow(sc)+1):(nrow(sc)*2)])
fits$residual_t = pars[length(pars)]
fits$halflife_t <- log(2)/fits$beta_t/60
fits$offset <- pars[length(pars)-1]
# print(fits)
hlo <- merge(hl,fits,by="accession")
print(cor(hlo$halflife,hlo$halflife_t,method="spearman"))
print(cor(hlo$alpha,hlo$alpha_t,method="spearman"))
write.table(fits,file=args[2],sep="\t",row.names=FALSE,quote=FALSE)
