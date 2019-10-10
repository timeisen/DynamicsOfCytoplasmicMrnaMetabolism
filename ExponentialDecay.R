#Read data in 
library(ggplot2)
library(tictoc)
args <- commandArgs(trailingOnly=TRUE)
# sc <- read.table("/lab/solexa_bartel/teisen/RNAseq/kinetics_of_translation/final_analyses/analysis/processed_files/miR-155_RNA_minus_tags_10rpm_cutoff_with_stds_scaled.txt",head=TRUE,sep="\t")
sc <- read.table(args[1],head=TRUE,sep="\t")
# sc <- head(sc,1000)
options(warn=1)

## Determine which dataset to use
if(args[3]=="miR-155_minus"){colnames(sc) <- 
	c("accession","x40m","x1h","x2h","x4h","x8h","xSS")
} else if(args[3]=="miR-155_plus"){colnames(sc) <- 
	c("accession","x40m","x1h","x2h","x4h","x8h","xSS") 
} else if(args[3]=="miR-1_minus"){colnames(sc) <- 
	c("accession","x40m","x1h","x2h","x4h","xSS") 
} else if(args[3]=="miR-1_plus"){colnames(sc) <-
	c("accession","x40m","x1h","x2h","x4h","xSS") 
} else {
	print("miRNA not specified")
	break} #miR-155 minus dox

## Remove the standards from the global fitting. 
stds <- read.table("/lab/solexa_bartel/teisen/RNAseq/Annotation_files/TAIL_seq/list_of_standards_SWE_notation.txt")
sc <- sc[which(!sc$accession %in% stds$V1),]

## Comparison to another dataset, to be used as a final sanity check.
hl<-read.table("/lab/solexa_bartel/eichhorn/5EU_RPF/miR-1_miR-155_analysis/Halflife_analysis/Halflife_measurements/miR-155_RNA_minusdox_halflife_SS.txt",head=TRUE,sep="\t")

## Main objective function
objective<-function(params,time,data,gene_no){
	alpha = exp(params[1:gene_no])
	beta = exp(params[(gene_no + 1):(gene_no*2)])
	offset = params[gene_no*2 + 1]

	scale = exp(params[gene_no*2 + 2])
	## Matrix to scale the final time point
	scalemat = matrix(c(rep(1,length(time) - 1),scale), 
		byrow = TRUE,nrow = gene_no,ncol = length(time))

	## Run the model
    model = sapply(time,function(t,alpha,beta,offset){
    				alpha/beta*(1-exp(-beta*(t-offset)))
					},
				alpha = alpha, 
				beta = beta, 
				offset=offset
				) * scalemat
	
    ## L2
	sqs <- sum((log(model)-log(data))^2)
	return(sqs)
}

## Analytic gradient for optimization efficiency

grAnalytic <- function(params,time,data,gene_no){
	alpha = exp(params[1:gene_no]) #because these are exponentiated here, I need to multiply the final "alphas" by alpha, same with beta.
	beta = exp(params[(gene_no+1):(gene_no*2)])
	offset = params[gene_no*2+1]
	
	scale = exp(params[gene_no*2 + 2])
	scalemat = matrix(c(rep(1,length(time) - 1),scale), 
		byrow = TRUE,nrow = gene_no,ncol = length(time))
	dscalemat = matrix(c(rep(0,length(time) - 1),1),
		byrow = TRUE,nrow = gene_no,ncol = length(time))

    model = sapply(time,function(t,alpha,beta,offset){
    	alpha/beta*(1-exp(-beta*(t-offset)))
    	},
    	alpha = alpha,
    	beta = beta,
    	offset = offset) * scalemat

	preFactor = 2 * (log(model) - log(data)) / model ## Derivative of loss

	## Derivative of alphas
	alphas = rowSums(preFactor * sapply(time,function(t,alpha,beta,offset){
		(1-exp(-beta*(t-offset)))/beta
		},
		alpha = alpha,
		beta = beta,
		offset = offset) * scalemat
		) * alpha
	
	## Derivative of betas
	betas = rowSums(preFactor * sapply(time,function(t,alpha,beta,offset){
		-alpha * (1-exp(-beta*(t-offset))) / beta^2 - 
		alpha * exp(-beta*(t-offset)) * (offset-t) / beta
		},
		alpha = alpha,
		beta = beta,
		offset = offset) * scalemat) * beta
	
	## Derivative of offset
	doffset = sum(preFactor * sapply(time,function(t,alpha,beta,offset){
		-alpha*exp(-beta*(t-offset))
		},
		alpha = alpha,
		beta = beta,
		offset = offset) * scalemat)

	## Derivative of the scale
	dscale = sum(preFactor*sapply(time,function(t,alpha,beta,offset){
		alpha / beta * (1-exp(-beta*(t-offset)))
		},
		alpha=alpha,
		beta=beta,
		offset=offset) * dscalemat
		) * scale

	## return the gradient
	gradient = c(alphas,betas,doffset,dscale)
    return(gradient)
}

## Runningn L-BFGS-B for optimization
optimize<-function(time,data,guess){
	gene_no = nrow(data)
	test_param <- c(log(rep(c(guess,1),each = gene_no)),32,log(1))
	for(i in 1:2){
	solve<-optim(
		par=test_param,
		fn=objective,
		gr=grAnalytic,
		method="L-BFGS-B",
		lower = c(log(rep(c(1E-10,1E-5),each = gene_no)),1,log(1E-1)),
		upper = c(log(rep(c(1E2,1E2),each = gene_no)),39,log(1E1)),
		data=data[-c(1)],
		time=time,
		gene_no = gene_no)
	test_param = solve$par}
	return(
		c(exp(solve$par[-length(solve$par)]),
		solve$par[length(solve$par) - 1], 
		exp(solve$par[length(solve$par)]), solve$value)
	)
}

if(args[3]=="miR-155_minus"){time_pts<-c(40,60,120,240,480,6000) #miR-155 minus dox
} else if(args[3]=="miR-155_plus"){
	time_pts <- c(41.0, 62.0, 120.0, 240.0, 480.0,6000) #miR-155 plus dox
} else if(args[3]=="miR-1_minus"){
	time_pts <- c(40.0, 70.0, 120.0, 241.0,6000) #miR-155 minus dox
} else if(args[3]=="miR-1_plus"){
	time_pts <- c(40.0, 65.0, 121.0, 240.0,6000) #miR-155 minus dox
} else {
	print("miRNA not specified")
	break} #miR-155 minus dox

guess = median(sc[,2])
tic("Optimize")
pars <- optimize(time_pts,sc,guess)
toc()

fits <- data.frame(
	accession = sc[,1],
	alpha_t = pars[1:nrow(sc)],
	beta_t = pars[(nrow(sc)+1):(nrow(sc)*2)]
	)

fits$residual_t = pars[length(pars)]
fits$halflife_t <- log(2)/fits$beta_t/60
## Metrics
print(median(fits$halflife_t))
print(median(fits$residual_t))

## Add in global params. 
fits$offset <- pars[length(pars)-2]
fits$scale <- pars[length(pars)-1]

## Comparison with previous data
hlo <- merge(hl,fits,by="accession")
print(cor(hlo$halflife,hlo$halflife_t,method="spearman"))
print(cor(hlo$alpha,hlo$alpha_t,method="spearman"))

## writing. 
write.table(fits,file=args[2],sep="\t",row.names=FALSE,quote=FALSE)

