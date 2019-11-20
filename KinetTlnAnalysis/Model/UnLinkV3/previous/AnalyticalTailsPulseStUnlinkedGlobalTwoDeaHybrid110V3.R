#populate vector of tails

args<-commandArgs(trailingOnly=TRUE)
library(deSolve)
library(numDeriv)
library(tictoc)
library(parallel)
library(data.table)
library(nloptr)
library(optimParallel)

# system("R CMD SHLIB /lab/solexa_bartel/teisen/RNAseq/Scripts/models/compiled/GlobalODE_10_plogis_linked.c")
# system("R CMD SHLIB /lab/solexa_bartel/teisen/RNAseq/Scripts/models/compiled/R_helpers.c")


# system("gcc -std=gnu99 -I/usr/share/R/include -DNDEBUG      -fpic -O3 -fstack-protector --param=ssp-buffer-size=4 -Wformat -Werror=format-security -D_FORTIFY_SOURCE=2 -c /lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/KinetTlnAnalysis/Model/UnLinkV2/AnalyticalTailsPulseStUnlinkedGlobalTwoDea110.c -o /lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/KinetTlnAnalysis/Model/UnLinkV2/AnalyticalTailsPulseStUnlinkedGlobalTwoDea110.o")
# system("gcc -std=gnu99 -shared -L/usr/lib/R/lib -Wl,-Bsymbolic-functions -Wl,-z,relro -o /lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/KinetTlnAnalysis/Model/UnLinkV2/AnalyticalTailsPulseStUnlinkedGlobalTwoDea110.so /lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/KinetTlnAnalysis/Model/UnLinkV2/AnalyticalTailsPulseStUnlinkedGlobalTwoDea110.o -L/usr/lib/R/lib -lR")

#mDCP has multiple deadenylation rates specified by a logistic function

dyn.load("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/KinetTlnAnalysis/Model/UnLinkV2/AnalyticalTailsPulseStUnlinkedGlobalTwoDea110.so") #load the model dynamically
dyn.load("/lab/solexa_bartel/teisen/RNAseq/Scripts/models/compiled/R_helpers_V2.so") #load the model dynamically

options(warn = 1)

#V9 allows all 4 parameters but incorporates expression into the measurements in order to determine a,b
#V9b fixes bugs in the list management
#V10b adopts this script to consider different optimal starting tail lengths.
#V11 allows only 3 parameters (a, k, b) but performs an optimization on a range of starting tail lengths in 10nt increments from 13 to 25. 
#V16 is a reversion to V11 that now tries to fit steady state tail length and increases the range that decapping can occur at to 90nt. 
#V20 uses matrices to solve the differential equations.
#V24 uses gradient optimization using the L-BFGS-B method, bounded constraints, and numDeriv gradient calculations. In addition, it implements a new system for measuring starting tail length using a single matrix and a starting tail length distribution determined by a gaussian with a sd of 1 and mean=stl. 
#the array script changes first lines of this file to make it compatible with a job array.

#initial param
#8nt trim from smallest tail lengths
#way up on transcription rate?
#logis
#V48 fits only one b scaling term
#V49 change the starting distribution to negative binomial
#57 has box constraints
#V59 uses randomized gene lists
#V64 uses a smoothing function and fits the global params c2, shape, size, etc. 
#V68 uses a logistic function for decapping

#read data in
#uncomment below to run
all_data_fn<-"/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/model_input_files/randomized_gene_sets/miR-155_minus_rand_g_v6_last8rm_HYBRID20190326_"
i <- args[1]
all_data<-fread(paste0(all_data_fn,i,".txt"),head=TRUE)
data<-as.matrix(all_data)

#uncomment below to test
# all_data<-fread(args[1],head=TRUE)
# all_data<-fread("model_input_files/TEST5_genes_miR-155_minus_sample_means_background_subtracted_v6_st_global_last8rm_HYBRID.txt",head=TRUE)
# data<-as.matrix(all_data)
# i<-"test"
# print(detectCores())
#dataset params
offset=35
time_points<-(c(40,60,120,240,480,6000)-offset)
#initial_param<-c(1.403647e+02,100,6.551659e-02,2.718854e-02) #what should these be?
N <<- ncol(data)/250 #number of genes
accession<<-colnames(all_data)[1+250*c(0:(N-1))]
initial_param<-c(rep(c(140,1E-7,1,1,1),N),7.2,241,10)
out_file <- paste0("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/rate_constant_measurements/global_model_fits/v7/Multicore_20190326_HYBRID_minus_TwoDea_110_",i,".txt")

data_split <- split((1:ncol(data)), ceiling((1:ncol(data))/250))
data <- unlist(lapply(data_split,function(index,data){data[,index]},data))

# minus
vars<<-rep(c(
  3.421246e-16,
  7.775342e-16,
  6.607287e-15,
  2.794067e-14,
  6.053765e-14,
  5.036801e-13),250*N)

# #plus
# vars<<-rep(c(
#   2.761815e-16,
#   1.310231e-15,
#   1.215492e-14,
#   2.101872e-14,
#   9.565653e-14,
#   9.550848e-13),250*N)



vars <- 1/vars #to allow multiplication in the c code.
###############################################################
##removes fitting values for non-steady state conditions.
oneGeneRm = c((6*242+1):(6*250))
oneGeneRm = oneGeneRm[-c((1:8)*6)]
val0 <- unlist(lapply(oneGeneRm,function(x){x+(0:(N-1))*6*250}))
vars[val0] = 0
###############################################################



Simulation<-function(pars,time){ ##cleanup
  max_tail = 251
  initial_state <- rep(0,max_tail)
  tails <- ode.band(func = "ode_deriv", y = initial_state, parms = pars, times=c(0,time),method = "lsode",bandup = 0, banddown = 1, nspec = max_tail, dllname = "AnalyticalTailsPulseStUnlinkedGlobalTwoDea110",nout=1,initfunc = "ode_p_init",hmax=60)
  columns_to_remove = c(1,2,ncol(tails))
  sim <- tails[-1,-columns_to_remove]
  return(sim)
}



tick <<- 0
CalculateResidual <- function(pars,data,time_points) {
  tick<<-tick+1
  param_list <- split(pars[1:(N*5)], ceiling(1:(N*5)/5))  
  # tic()
  # no_cores <- detectCores()
  # cl <- makeCluster(no_cores)
  model <- unlist(mclapply(param_list,function(gene_par,pars,time_points){Simulation(2^c(gene_par,tail(pars,3)),time_points)},pars,time_points))
  residual <- .C("ReturnResidualMod",n = as.integer(N),model = as.double(model),data = as.double(data),var = as.double(vars),residual = as.double(0))$residual
  # stopCluster(cl)
  # toc()

  # if(tick%%10==0){
  #   print(tick/10)
  #   print(residual)}
  if(is.finite(residual)){return(residual)}
  else(return(10E20))
  return(residual)
}

grr<-function(pars,data,time_points){
    gradient<-grad(CalculateResidual,pars,data=data,method="simple",time_points=time_points)
    return(gradient)
}


Optimization<-function(initial_param,data,plot=FALSE,time_points){
  #print(paste("starting_tail_length",starting_tail_length))
  print(length(initial_param))
  if(plot){plot(time_points,log(data[6:10]), pch = 19,ylim=c(0,10),xlim=c(0,900))} #only plotting tails
  solve<-NULL
  optim_param<-NULL
  for(x in 1:1){
    solve$solution<-initial_param
    for(i in 1:3){
      solve<-nloptr(
      x0=solve$solution,
      eval_f=CalculateResidual,
      # localsolver = c("lbfgs"),
      # method="L-BFGS-B",
      eval_grad_f = grr,
      data=data,
      # plot=FALSE,
      lb=log2(c(rep(c(10E-10,10E-10,10E-10,10E-10,10E-10),N),0,0,0)),
      ub=log2(c(rep(c(10E10,10E10,10E10,10E10,10E10),N),100,500,500)),
      time_points=time_points,
      # nl.info=TRUE,
      opts=list("algorithm"="NLOPT_LD_LBFGS",xtol_rel=1e-8,maxeval=300,"print_level"=1))
      print(solve$objective)
      }
    optim_param<-rbind(optim_param,c(2^solve$solution,solve$objective))
    }
  return(optim_param)
  }


all_optimizations<-NULL
tic("optim_c")
all_optimizations<-Optimization(log2(initial_param),data,plot=FALSE,time=time_points)
toc()

lens <- length(all_optimizations)
optim_data<-data.frame(matrix(all_optimizations[-c(lens,lens-1,lens-2,lens-3)],ncol=5,byrow=TRUE))
colnames(optim_data)<-c("st","a","k1","k2","b")
optim_data$size <- all_optimizations[lens-3]
optim_data$location <- all_optimizations[lens-2]
optim_data$scale <- all_optimizations[lens-1]
optim_data$residual <- all_optimizations[lens]
optim_data$accession <- accession

# print(optim_data)
print(tick)
write.table(optim_data,file=out_file,row.names=FALSE,quote=FALSE,sep="\t",append=FALSE)
# write.table(final_set,file="test.txt",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t",append=FALSE)
