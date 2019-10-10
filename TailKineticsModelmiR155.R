args<-commandArgs(trailingOnly=TRUE)
options(warn = 1)

## packages and imports
library(deSolve)
library(numDeriv)

## Main model is written in c and loaded dynamically.
## It can be compiled with the next command. 
#system("R CMD SHLIB AnalyticalTailsPulseStUnlinked.c")
dyn.load("AnalyticalTailsPulseStUnlinked.so") 

all_data <- read.table(args[1],head=TRUE)
accession <<- as.character(all_data[args[2],1])
data <- as.numeric(all_data[args[2],-1])
offset = 35 #Minutes to offset the data to account for export?
time_points <- (c(40,60,120,240,480,6000)-offset)
initial_param <- c(140,1E-7,1,1)
dir <- args[3]

##Average the last 20 nt for fitting##
data[(length(data) - 19):length(data)] = 
  mean(data[(length(data) - 19):length(data)])
##

#The variances of the datasets, to be used for residual weighting. 
#Global parameter assignment. 

vars <<- c(      #miR-155 minus
  7.426524e-16,
  1.287242e-15,
  3.640165e-15,
  9.202623e-15,
  4.320013e-14,
  2.808383e-13,
  2.808383e-13/6) 
## This increases the weighting value of these points 6 fold
## Updated variances as of 2019 06 03.


Simulation <- function(pars,time){ ##This is the main simulator

  #Parameter definitions
  st       = pars[1]
  a        = pars[2]
  k        = pars[3]
  b        = pars[4]
  size     = 16.25643 #Global fitted parameters as of 2019 07 15.
  location = 263.95156 
  scale    = 11.05133 

  parameters = c(st,a,k,b,size,location,scale)
  max_tail = 251

  initial_state <- rep(0,max_tail) #All abundances begin with 0. 
  #The simulation, passed to c code called ode_deriv, using lsode.
  #This is for a banded jacobian.
  #Note hmax has a major impact on memory usage, time, and precision. 
  tails <- ode.band(func = "ode_deriv", y = initial_state, parms = parameters, 
          times = c(0,time), method = "lsode",bandup = 0, banddown = 1,
          nspec = max_tail, dllname = "AnalyticalTailsPulseStUnlinked",
          nout=1, initfunc = "ode_p_init", hmax = 1,maxsteps=5000000)

  ##These two lines below return NA if the tails output is incomplete. Important
  # for using randomized initial parameters.
  if(dim(tails)[1] == 7){}
  else(return(NA))
  
  #Remove columns that shouldn't be compared to residuals. 
  #Why is it 18? ncol(tails) = 253, max_tail = 251, 233:253 is 
  #      21 values including the last NA column
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

  #Dealing with errors in residual calculation.
  if(is.finite(residual)){return(residual)}
  else(return(10E20))
}

#Simple, finite differences gradient calculation.
grr <- function(pars,data,plot=plot,time_points){
    gradient <- grad(CalculateResidual, pars, 
      data=data, method="simple",time_points=time_points)
    return(gradient)
}

#Using L-BFGS-B algorithm
Optimization <- function(initial_param, data, plot=FALSE, time_points){
  #print(paste("starting_tail_length",starting_tail_length))
  if(plot){
        plot(time_points,log(data[6:10]), 
        pch = 19,ylim=c(0,10),xlim=c(0,900))} #only plotting tails
  solve <- NULL
  optim_param <- NULL
  for(x in 1:2){
    initial_param <- log2(runif(4,
      min=c(147.16777,4.317860e-08,2.650033e-01,2.650033e-01),
      max=c(168.41418,1.396230e-07,9.549402e-01,9.549402e-01)))
    # scale<-CalculateParscale(initial_param,time_points)
    solve$par<-initial_param
    for(i in 1:2){
      #arguments passed to the solver
      solve<-optim(
        par=solve$par,
        fn=CalculateResidual,
        method="L-BFGS-B",
        gr = grr,
        data=data,
        lower=log2(c(1,0,0,0)),
        upper=log2(c(250,Inf,Inf,Inf)),
        time_points=time_points,
        plot=plot,
        control=list())}
    optim_param<-rbind(optim_param,c(2^solve$par,solve$value))
    }
  return(optim_param)
  }

#Main code block for running the optimimization and writing output files. 
all_optimizations <- NULL
all_optimizations <- Optimization(
    log2(initial_param),
    data,plot=FALSE, time=time_points)
# print(accession)
# print(all_optimizations)
all_optimizations<- cbind( 
  rep(accession,nrow(all_optimizations)),
  all_optimizations)
fn <- paste0("OUTPATH/",
  dir,"/",accession,".txt")
write.table(all_optimizations,
  file=fn,
  row.names=FALSE,
  col.names=FALSE,
  quote=FALSE,
  sep="\t",
  append=FALSE)

