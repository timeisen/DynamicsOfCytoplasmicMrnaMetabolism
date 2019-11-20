#populate vector of tails

args<-commandArgs(trailingOnly=TRUE)
library(deSolve)
library(numDeriv)
library(tictoc)
# system("R CMD SHLIB /lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/KinetTlnAnalysis/Model/AnalyticalTailsPulseStUnlinked.c")
dyn.load("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/KinetTlnAnalysis/Model/UnLinkV2/AnalyticalTailsPulseStUnlinkedGlobalTwoDea150.so") #load the model dynamically
options(warn = 1)
##model_input_files/miR-155_minus_sample_mean_tails_50tags_with_expr_background_subtracted_v7_st_last7rm_HYBRID20190326.txt

#V9 allows all 4 parameters but incorporates expression into the measurements in order to determine a,b
#V9b fixes bugs in the list management
#V10b adopts this script to consider different optimal starting tail lengths.
#V11 allows only 3 parameters (a, k, b) but performs an optimization on a range of starting tail lengths in 10nt increments from 13 to 25. 
#V16 is a reversion to V11 that now tries to fit steady state tail length and increases the range that decapping can occur at to 90nt. 
#V20 uses matrices to solve the differential equations.
#V24 uses gradient optimization using the L-BFGS-B method, bounded constraints, and numDeriv gradient calculations. In addition, it implements a new system for measuring starting tail length using a single matrix and a starting tail length distribution determined by a gaussian with a sd of 1 and mean=stl. 
#the array script changes first lines of this file to make it compatible with a job array.
#V75 is unlinked

#initial param
#8nt trim from smallest tail lengths
#way up on transcription rate?
#pexp
#V48 fits only one b scaling term
#V49 change the starting distribution to negative binomial
#57 has box constraints
#V58 with smoothing
#V71 uses a plogis function in the ode model. This is the current version of the script as of 20180221

#linked: 213502562
#unlinked: 225708248

all_data<-read.table(args[1],head=TRUE)
accession<<-as.character(all_data[args[2],1])
data<-as.numeric(all_data[args[2],-1])
offset=35
time_points<-(c(40,60,120,240,480,6000)-offset)
#initial_param<-c(1.403647e+02,100,6.551659e-02,2.718854e-02) #what should these be?
initial_param<-c(140,1E-7,1,1,1)
dir <- args[3]
N<<-1
print(accession)
##Average the last 20 nt for fitting##
data[(length(data) - 19):length(data)] = mean(data[(length(data) - 19):length(data)])
##

vars <<- c(      #miR-155 minus
  7.426524e-16,
  1.287242e-15,
  3.640165e-15,
  9.202623e-15,
  4.320013e-14,
  2.808383e-13,
  2.808383e-13/6) #The variance of the last value is reduced 6 fold
## This increases the weighting value of these points 6 fold
#Updated variances as of 2019 06 03.

Simulation<-function(pars,time){ ##cleanup

  st       = pars[1]
  a        = pars[2]
  k1       = pars[3]
  k2       = pars[4]
  b        = pars[5]
  size     = 12.8628040    #110 21.1387083, 25 datasets. #150: 13.8574646   
  location = 264.1838832   #110 260.0190246 #150: 260.4219318 
  scale    = 10.0829122    #110 9.1010557 #150: 8.7416984   

  parameters = c(st,a,k1,k2,b,size,location,scale)
  max_tail = 251
  initial_state <- rep(0,max_tail*N)

  tails <- ode.band(func = "ode_deriv", 
    y = initial_state, 
    parms = parameters, 
    times=c(0,time),
    method = "lsode",
    bandup = 0, 
    banddown = 1,
    nspec = max_tail*N, 
    dllname = "AnalyticalTailsPulseStUnlinkedGlobalTwoDea150",
    nout=1,
    initfunc = "ode_p_init",
    hmax = 1,
    maxsteps=5000000)

  #Remove columns that shouldn't be compared to residuals. 
  #Why is it 18? ncol(tails) = 253, max_tail = 251, 233:253 is 21 values including the last NA column
  #      which puts it at 20 values including the 0.
  columns_to_remove = c(1,2,(max_tail - 18):ncol(tails))
  #Note that these values are all replaced by the mean here for the last 20 nt.
  last20ntMean = rep(mean(tails[7,(max_tail - 18):(ncol(tails)-1)]), 20)
  sim <- tails[-1,-columns_to_remove]
  #Add the last 20 nt back to the flattened array. 
  sim <- c(c(t(sim)),last20ntMean)
  return(sim)
}

tick = 0
CalculateResidual <- function(pars,data,plot=FALSE,time_points) { #Run the model

  model<-Simulation(2^pars,time_points) #exponentiation the parameters.

  #Residuals, weighted. 
  sqs <- (model-data)^2
  sqs[1:230]     <- sqs[1:230]    / vars[1]
  sqs[231:460]   <- sqs[231:460]  / vars[2]
  sqs[461:690]   <- sqs[461:690]  / vars[3]
  sqs[691:920]   <- sqs[691:920]  / vars[4]
  sqs[921:1150]  <- sqs[921:1150] / vars[5]
  sqs[1151:1380] <- sqs[1151:1380]/ vars[6]
  sqs[1381:1400] <- sqs[1381:1400]/ vars[7] #Added in V4 by TJE on 2019 04 01. 

  residual <- sum(sqs*1E6)
  
  #For plotting
  if (tick%%10 == 0 & plot){
    final_set <<- matrix(c(model,data,rep(tick,1500)),ncol=3)
  }
  tick <<- tick + 1
  # print(residual)
  #Dealing with errors in residual calculation.
  if(is.finite(residual)){return(residual)}
  else(return(10E20))
}


grr<-function(pars,data,plot=plot,time_points){
    gradient<-grad(CalculateResidual,pars,data=data,method="simple",time_points=time_points)
    return(gradient)
}

Optimization<-function(initial_param,data,plot=FALSE,time_points){
  #print(paste("starting_tail_length",starting_tail_length))
  if(plot){plot(time_points,log(data[6:10]), pch = 19,ylim=c(0,10),xlim=c(0,900))} #only plotting tails
  solve<-NULL
  optim_param<-NULL
  for(x in 1:1){
    initial_param <- log2(runif(5, min=c(134.2714, 2.366310e-08, 0.1321663, 0.1321663, 0.03031561), max=c(188.2605,6.922078e-07,6.4182735,6.4182735,26.08015665))) #this is 5th to 95th percentile
    # scale<-CalculateParscale(initial_param,time_points)
    solve$par<-initial_param
    for(i in 1:2){
      solve<-optim(
      par=solve$par,
      fn=CalculateResidual,
      method="L-BFGS-B",
      gr = grr,
      data=data,
      lower=log2(c(1,0,0,0,0)),
      upper=log2(c(250,Inf,Inf,Inf,Inf)),
      time_points=time_points,
      plot=plot,
      control=list())}
    optim_param<-rbind(optim_param,c(2^solve$par,solve$value))
    }
  return(optim_param)
  }


all_optimizations<-NULL
all_optimizations<-Optimization(log2(initial_param),data,plot=FALSE,time=time_points)
# print(accession)
# print(all_optimizations)
all_optimizations<-cbind(rep(accession,nrow(all_optimizations)),all_optimizations)
fn<-paste("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/optim_runs/",dir,"/",accession,".txt",sep="")
write.table(all_optimizations,file=fn,row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t",append=FALSE)

