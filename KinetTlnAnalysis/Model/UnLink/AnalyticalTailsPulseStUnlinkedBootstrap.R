#populate vector of tails

args<-commandArgs(trailingOnly=TRUE)
library(deSolve)
library(numDeriv)
library(tictoc)
# system("R CMD SHLIB /lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/KinetTlnAnalysis/Model/AnalyticalTailsPulseStUnlinked.c")
dyn.load("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/KinetTlnAnalysis/Model/AnalyticalTailsPulseStUnlinked.so") #load the model dynamically
options(warn = 1)

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
bootstrap_params <- read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/rate_constant_measurements/global_model_fits/bootstrapUnlinked/miR-155_unlinked_rates_compiled_median_20180709.txt",sep="\t",head=TRUE)
global_params <<- as.numeric(bootstrap_params[args[4],2:4])

all_data<-read.table(args[1],head=TRUE)
accession<<-as.character(all_data[args[2],1])
data<-as.numeric(all_data[args[2],-1])
offset=35
time_points<-(c(40,60,120,240,480,6000)-offset)
#initial_param<-c(1.403647e+02,100,6.551659e-02,2.718854e-02) #what should these be?
initial_param<-c(140,1E-7,1,1)
dir <- args[3]
N<<-1

vars<<-c(      #miR-155 minus
  2.58757e-16,
  3.053586e-16,
  2.92142e-15,
  1.16722e-14,
  2.772866e-14,
  2.394499e-13)

# vars<<-c(      #miR-155 plus
#   NA1.57175e-16,
#   NA3.988673e-16,
#   NA2.758914e-15,
#   NA5.179771e-15,
#   NA1.93081e-14,
#   NA1.569503e-13)


Simulation<-function(pars,time){ ##cleanup

  st       = pars[1]
  a        = pars[2]
  k        = pars[3]
  b        = pars[4]
  size     = global_params[1]
  location = global_params[2]
  scale    = global_params[3]

  parameters = c(st,a,k,b,size,location,scale)
  max_tail = 251
  initial_state <- rep(0,max_tail*N)
  tails <- ode.band(func = "ode_deriv", y = initial_state, parms = parameters, times=c(0,time),method = "lsode",bandup = 0, banddown = 1,nspec = max_tail*N, dllname = "AnalyticalTailsPulseStUnlinked",nout=1,initfunc = "ode_p_init",hmax = 60)

  columns_to_remove = c(1,2,ncol(tails))
  sim <- tails[-1,-columns_to_remove]
  return(c(t(sim)))

}

tick = 0
CalculateResidual <- function(pars,data,plot=FALSE,time_points) {
  model<-Simulation(2^pars,time_points)
  sqs<-(model-data)^2
  sqs[1:250]<-sqs[1:250]/vars[1]
  sqs[251:500]<-sqs[251:500]/vars[2]
  sqs[501:750]<-sqs[501:750]/vars[3]
  sqs[751:1000]<-sqs[751:1000]/vars[4]
  sqs[1001:1250]<-sqs[1001:1250]/vars[5]
  sqs[1251:1500]<-sqs[1251:1500]/vars[6]
  residual <- sum(sqs*1E6)
  if (tick%%10 == 0 & plot){
    final_set <<- matrix(c(model,data,rep(tick,1500)),ncol=3)
  }
  tick <<- tick + 1
  if(is.finite(residual)){return(residual)}
  else(return(10E20))
}

CalculateParscale<-function(pars,time_points){
  scale<-pars
  data<-Simulation(2^pars,time_points)
  for(i in 1:length(pars)){
    newpars<-pars
    newpars[i]<-pars[i]+0.1
    model<-Simulation(2^(newpars),time_points)
    scale[i]<-sum((model-data)^2)
    }
  return(scale)
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
    initial_param<-log2(runif(4,min=c(147.16777,4.317860e-08,2.650033e-01,100*2.650033e-01),max=c(168.41418,1.396230e-07,9.549402e-01,100*9.549402e-01)))
    # scale<-CalculateParscale(initial_param,time_points)
    solve$par<-initial_param
    for(i in 1:2){
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


all_optimizations<-NULL
all_optimizations<-Optimization(log2(initial_param),data,plot=FALSE,time=time_points)
# print(accession)
# print(all_optimizations)
all_optimizations<-cbind(rep(accession,nrow(all_optimizations)),all_optimizations)
fn<-paste("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/optim_runs/",dir,"/",accession,".txt",sep="")
write.table(all_optimizations,file=fn,row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t",append=FALSE)

