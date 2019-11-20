#populate vector of tails

args<-commandArgs(trailingOnly=TRUE)
library(deSolve)
library(numDeriv)
library(tictoc)
library(parallel)

########## DOES NOT RUN YET #############
########## Warning in lsode(y, times, func, parms = parms, bandup = bandup, banddown = banddown,  :
##########   repeated convergence test failures on a step, but integration was successful - inaccurate Jacobian matrix?
########## Warning in lsode(y, times, func, parms = parms, bandup = bandup, banddown = banddown,  :
##########   Returning early. Results are accurate, as far as they go
########## Warning in mclapply(data, Optimization, initial_param = initial_param, time_points = time_points,  :
##########   scheduled core 1 encountered error in user code, all values of the job will be affected
########## Testing: 10194.617 sec elapsed
########## Warning in matrix(unlist(all_optimizations), byrow = TRUE, ncol = 5) :
##########   data length [9247] is not a sub-multiple or multiple of the number of rows [1850]
########## Error in data.frame(..., check.names = FALSE) : 
##########   arguments imply differing number of rows: 3083, 1850
########## Calls: cbind -> cbind -> data.frame
########## Execution halted


# system("R CMD SHLIB /lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/KinetTlnAnalysis/Model/AnalyticalTailsPulseStUnlinked.c")
dyn.load("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/KinetTlnAnalysis/Model/UnLink/AnalyticalTailsPulseStUnlinked.so") #load the model dynamically
dyn.load("/lab/solexa_bartel/teisen/RNAseq/Scripts/models/compiled/R_helpers_V3.so") #residual sum scripts
options(warn = 1)
print(paste(c("cores:",detectCores())),quote=FALSE)

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
data <- split(all_data[,-1], seq(nrow(all_data[-1])))
offset=35
time_points<-(c(40,60,120,240,480,6000)-offset)
#initial_param<-c(1.403647e+02,100,6.551659e-02,2.718854e-02) #what should these be?
initial_param<-c(140,1E-7,1,1)
N<-1

vars<-rep(c( ##miR-155_minus
  2.66568e-16,
  3.142525e-16,
  3.00228e-15,
  1.200395e-14,
  2.854782e-14,
  2.470686e-13),each = 242)


vars[242*(1:(N*6))] <- vars[242*(1:(N*6))]*8 #adjusts the variance for those 8 nt. 
vars <- 1/vars #to allow multiplication in the c code.
print(length(vars))



Simulation<-function(pars,time){ ##cleanup

  st       = pars[1]
  a        = pars[2]
  k        = pars[3]
  b        = pars[4]
  size     = 12
  location = 286
  scale    = 12.7
  N = 1
  parameters = c(st,a,k,b,size,location,scale)
  max_tail = 251
  initial_state <- rep(0,max_tail*N)
  tails <- ode.band(func = "ode_deriv", y = initial_state, parms = parameters, times=c(0,time),method = "lsode",bandup = 0, banddown = 1,nspec = max_tail*N, dllname = "AnalyticalTailsPulseStUnlinked",nout=1,initfunc = "ode_p_init",hmax = 60)

  columns_to_remove = c(1,2,ncol(tails))
  sim <- tails[-1,-columns_to_remove]
  return(sim)

}

CalculateResidual <- function(pars,data,plot=FALSE,time_points,vars) {
  sim<-Simulation(2^pars,time_points)
  if(is.matrix(sim)){
    model <- c(t(cbind(sim[,1:241],rowSums(sim[,242:250])))) ##This line sums the last 8 nt. 
    residual <- .C("ReturnResidualModScale",n = as.integer(N),model = as.double(model),data = as.double(data),var = as.double(vars),residual = as.double(0))$residual*1E6
    return(residual)}
  else(return(10E20))
}

grr<-function(pars,data,plot=plot,time_points,vars){
    gradient<-grad(CalculateResidual,pars,data=data,method="simple",time_points=time_points,vars=vars)
    return(gradient)
}

Optimization<-function(data,initial_param,plot=FALSE,time_points,vars){
  #print(paste("starting_tail_length",starting_tail_length))
  if(plot){plot(time_points,log(data[6:10]), pch = 19,ylim=c(0,10),xlim=c(0,900))} #only plotting tails
  solve<-NULL
  optim_param<-NULL
  for(x in 1:1){
    initial_param<-log2(runif(4,min=c(147.16777,4.317860e-08,2.650033e-30,100*2.650033e-01),max=c(168.41418,1.396230e-07,9.549402e-01,100*9.549402e-01)))
    solve$par<-initial_param
    for(i in 1:2){
      solve<-optim(
      par=solve$par,
      fn=CalculateResidual,
      method="L-BFGS-B",
      gr = grr,
      data=data,
      vars = vars,
      lower=log2(c(1,0,0,0)),
      upper=log2(c(250,Inf,Inf,Inf)),
      time_points=time_points,
      plot=plot,
      control=list())}
    optim_param<-rbind(optim_param,c(2^solve$par,solve$value))
    }
  return(optim_param)
  }


all_optimizations<-NULL
tic("Testing")
all_optimizations<-mclapply(data,Optimization,initial_param = initial_param,time_points = time_points, vars = vars)
toc()
all_optimizations_df <- data.frame(matrix(unlist(all_optimizations),byrow=TRUE,ncol=5))
all_optimizations_df <- cbind(all_data$accession,all_optimizations_df)

write.table(all_optimizations_df,file=args[2],row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t",append=FALSE)

