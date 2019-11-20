#Read data in 
library(ggplot2)
library(tictoc)
args <- commandArgs(trailingOnly=TRUE)
# sc <- read.table("/lab/solexa_bartel/teisen/RNAseq/kinetics_of_translation/final_analyses/analysis/processed_files/miR-155_RNA_minus_tags_10rpm_cutoff_with_stds_scaled.txt",head=TRUE,sep="\t")
sc <- read.table(args[1],head=TRUE,sep="\t")
options(warn=1)

# sc <- sc[,c(1,8:13)]
# sc <- head(sc,500)
if(args[3]=="miR-155_minus"){colnames(sc)<-c("accession","x40m","x1h","x2h","x4h","x8h","xSS")
} else if(args[3]=="miR-155_plus"){colnames(sc)<-c("accession","x40m","x1h","x2h","x4h","x8h","xSS") 
} else if(args[3]=="miR-1_minus"){colnames(sc)<-c("accession","x40m","x1h","x2h","x4h","xSS") 
} else if(args[3]=="miR-1_plus"){colnames(sc)<-c("accession","x40m","x1h","x2h","x4h","xSS") 
} else {
	print("miRNA not specified")
	break} #miR-155 minus dox

stds<-read.table("/lab/solexa_bartel/teisen/RNAseq/Annotation_files/TAIL_seq/list_of_standards_SWE_notation.txt")

sc <- sc[which(!sc$accession %in% stds$V1),]

hl<-read.table("/lab/solexa_bartel/eichhorn/5EU_RPF/miR-1_miR-155_analysis/Halflife_analysis/Halflife_measurements/miR-155_RNA_minusdox_halflife_SS.txt",head=TRUE,sep="\t")
incorObjective = 0
objective<-function(params,time,data,gene_no){
	k0 = exp(params[1:gene_no])
	k1 = exp(params[gene_no*2+1])
	k2 = exp(params[(gene_no+1):(gene_no*2)])
    model = sapply(time,function(t,k0,k1,k2){k0/k2+k0*exp(-k1*t)/(k1-k2)-k0*k1*exp(-k2*t)/((k1-k2)*k2)},k0=k0,k1=k1,k2=k2)
	sqs<-sum((log(model)-log(data))^2)
	if(!is.finite(sqs)){
		incorObjective = incorObjective+1
		return(1E8)}
	return(sqs)
}
incorGrad = 0
###Analytic gradient for optimization efficiency
grAnalytic <- function(params,time,data,gene_no){
	k0 = exp(params[1:gene_no]) #because these are exponentiated here, I need to multiply the final "alphas" by alpha, same with beta.
	k2 = exp(params[(gene_no+1):(gene_no*2)])
	k1 = exp(params[gene_no*2+1])
    model = sapply(time,function(t,k0,k1,k2){k0/k2+k0*exp(-k1*t)/(k1-k2)-k0*k1*exp(-k2*t)/((k1-k2)*k2)},k0=k0,k1=k1,k2=k2)
	preFactor = 2*(log(model)-log(data))/model
	#derivative from alphas
	k0gr = rowSums(preFactor*sapply(time,function(t,k0,k1,k2){exp(-k1*t)/(k1-k2)+1/k2-exp(-k2*t)*k1/((k1-k2)*k2)},k0=k0,k1=k1,k2=k2))*k0
	#derivative of offset
	k1gr = sum(preFactor*sapply(time,function(t,k0,k1,k2){-((exp(-k1*t)*k0)/(k1 - k2)^2) + (exp(-k2*t)*k0*k1)/((k1 - k2)^2*k2) - (exp(-k2*t)*k0)/((k1-k2)*k2) - (exp(-k1*t)*k0*t)/(k1 - k2)},k0=k0,k1=k1,k2=k2))*k1
	#derivative of betas
	k2gr = rowSums(preFactor*sapply(time,function(t,k0,k1,k2){(exp(-k1*t)*k0)/(k1 - k2)^2 - k0/k2^2 + (exp(-k2*t)*k0*k1)/((k1 - k2)*k2^2) - (exp(-k2*t)*k0*k1)/((k1-k2)^2*k2)+(exp(-k2*t)*k0*k1*t)/((k1 - k2)*k2)},k0=k0,k1=k1,k2=k2))*k2

	gradient = c(k0gr,k2gr,k1gr)
	if(any(!is.finite(gradient))){
		incorGrad = incorGrad+1
		return(runif(gene_no*2+1))} #this is a strange line
    return(gradient)
}

optimize<-function(time,data,guess){
	gene_no = nrow(data)
	test_param <- log(c(rep(c(guess,.02),each = gene_no),1))
	for(i in 1:2){
	solve<-optim(
		par=test_param,
		fn=objective,
		gr=grAnalytic,
		method="L-BFGS-B",
		lower = log(c(rep(c(1E-10,1E-5),each = gene_no),1E-5)),
		upper = log(c(rep(c(1E4,1E4),each = gene_no),1E4)),
		data=data[-c(1)],
		time=time,
		gene_no = gene_no)
	test_param = solve$par}
	return(c(exp(solve$par),solve$value))}

if(args[3]=="miR-155_minus"){time_pts<-c(40,60,120,240,480,6000) #miR-155 minus dox
} else if(args[3]=="miR-155_plus"){time_pts<-c(41.0, 62.0, 120.0, 240.0, 480.0,6000) #miR-155 plus dox
} else if(args[3]=="miR-1_minus"){time_pts<-c(40.0, 70.0, 120.0, 241.0,6000) #miR-155 minus dox
} else if(args[3]=="miR-1_plus"){time_pts<-c(40.0, 65.0, 121.0, 240.0,6000) #miR-155 minus dox
} else {
	print("miRNA not specified")
	break} #miR-155 minus dox
guess = median(sc[,2])
tic("Optimize")
pars <- optimize(time_pts,sc,guess)
toc()
print("Incorrect Objective:")
print(incorGrad)
print("Incorrect Gradient:")
print(incorObjective)
fits <- data.frame(accession = sc[,1],alpha_t = pars[1:nrow(sc)],beta_t = pars[(nrow(sc)+1):(nrow(sc)*2)])
fits$residual_t = pars[length(pars)]
fits$halflife_t <- log(2)/fits$beta_t/60
fits$k1 <- pars[length(pars)-1]
# print(fits)
hlo <- merge(hl,fits,by="accession")
print(cor(hlo$halflife,hlo$halflife_t,method="spearman"))
print(cor(hlo$alpha,hlo$alpha_t,method="spearman"))
write.table(fits,file=args[2],sep="\t",row.names=FALSE,quote=FALSE)
