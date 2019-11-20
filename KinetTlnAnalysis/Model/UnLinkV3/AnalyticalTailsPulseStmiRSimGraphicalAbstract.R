#populate vector of tails

library(deSolve)
library(numDeriv)
library(tidyverse)
# system("R CMD SHLIB /lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/KinetTlnAnalysis/Model/UnLinkV3/AnalyticalTailsPulseStmiRSim.c")
dyn.load("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/KinetTlnAnalysis/Model/UnLink/AnalyticalTailsPulseStUnlinked.so")
options(warn = 1)
source("/lab/solexa_bartel/teisen/RNAseq/Scripts/general/ggplot_theme.R")

# V9 allows all 4 parameters but incorporates expression into the measurements in order to determine a,b
# V9b fixes bugs in the list management
# V10b adopts this script to consider different optimal starting tail lengths.
# V11 allows only 3 parameters (a, k, b) but performs an optimization on a range of starting tail lengths in 10nt increments from 13 to 25. 
# V16 is a reversion to V11 that now tries to fit steady state tail length and increases the range that decapping can occur at to 90nt. 
# V20 uses matrices to solve the differential equations.
# V24 uses gradient optimization using the L-BFGS-B method, bounded constraints, and numDeriv gradient calculations. In addition, it implements a new system for measuring starting tail length using a single matrix and a starting tail length distribution determined by a gaussian with a sd of 1 and mean=stl. 
# the array script changes first lines of this file to make it compatible with a job array.
# V75 is unlinked
# 2018 09 18: This version removes the last 8 nt from the fitting and analysis.
# initial param
# 8nt trim from smallest tail lengths
# way up on transcription rate?
# pexp
# V48 fits only one b scaling term
# V49 change the starting distribution to negative binomial
# 57 has box constraints
# V58 with smoothing
# V71 uses a plogis function in the ode model. This is the current version of the script as of 20180221
# 2018 09 24 This version of the model uses global rate constants from the datasets that has 8 nt removed. It uses global parameters from a fitting that includes TAIL-seq data for steady state.
# The HYBRID code uses the last 8 nt from the TAIL-seq dataset to fit 850 genes.  
# 2019 03 01 The V3 is exactly the same as the V2 code, but run with new global parameters, varainces, and using the scaled values from the new data. 
# The V4 script files are updated for increasing the last 8 nt weighting 6 fold to account for the fact that we don't have those values in the steady state. 


#The variances of the datasets, to be used for residual weighting. 
#Global parameter assignment. 
#The variance of the last value is reduced 6 fold
## This increases the weighting value of these points 6 fold


cdff <- function(vec){
  vec <- vec[length(vec):1] #flip so it's 0:250
  vecCDF <- unlist(lapply(1:length(vec),function(x){sum(vec[1:x])}))
  vecCDF <- vecCDF/vecCDF[length(vec)] #normalize it
  return(vecCDF)
}

Simulation <- function(pars,time){ ##This is the main simulator

  #Parameter definitions
  st             = pars[1]
  a              = pars[2]
  k              = pars[3]
  b              = pars[4]
  size           = pars[5] #These params are from the prelim run, 2019 03 02, 22 datasets. 
  location       = pars[6] 
  scale          = pars[7] 

  parameters = c(st,a,k,b,size,location,scale)
  max_tail = 251

  initial_state <- rep(0,max_tail) #All abundances begin with 0. 

  #The simulation, passed to c code called ode_deriv, using lsode.
  #This is for a banded jacobian.
  #Note hmax has a major impact on memory usage, time, and precision. 
  tails <- ode.band(func = "ode_deriv", y = initial_state, parms = parameters, 
          times = c(0,time), method = "lsode",bandup = 0, banddown = 1,
          nspec = max_tail, dllname = "AnalyticalTailsPulseStUnlinked",
          nout=1, initfunc = "ode_p_init", hmax = 1, maxsteps=5000000)

  ##These two lines below return NA if the tails output is incomplete. Important
  # for using randomized initial parameters.
  
  #Remove columns that shouldn't be compared to residuals. 
  columns_to_remove = c(ncol(tails))
  sim <- tails[-1,-columns_to_remove]
  return(sim)

}
 

SimulationActD <- function(pars,time, is){ ##This is the main simulator

  #Parameter definitions
  st             = pars[1]
  a              = pars[2]
  k              = pars[3]
  b              = pars[4]
  size           = pars[5] #These params are from the prelim run, 2019 03 02, 22 datasets. 
  location       = pars[6] 
  scale          = pars[7] 

  parameters = c(st,a,k,b,size,location,scale)
  max_tail = 251

  initial_state <- is #All abundances begin with 0. 

  #The simulation, passed to c code called ode_deriv, using lsode.
  #This is for a banded jacobian.
  #Note hmax has a major impact on memory usage, time, and precision. 
  tails <- ode.band(func = "ode_deriv", y = is, parms = parameters, 
          times = c(0,time), method = "lsode",bandup = 0, banddown = 1,
          nspec = max_tail, dllname = "AnalyticalTailsPulseStUnlinked",
          nout=1, initfunc = "ode_p_init", hmax = 1, maxsteps=5000000)

  ##These two lines below return NA if the tails output is incomplete. Important
  # for using randomized initial parameters.
  
  #Remove columns that shouldn't be compared to residuals. 
  columns_to_remove = c(ncol(tails))
  sim <- tails[-1,-columns_to_remove]
  return(sim)

}

miR1SimPars1 <- c(
  180,
  1E-7,
  0.5/10,
  0.1/10,
  15.5,
  260,
  10)

#plus miR
miR1SimPars2 <- miR1SimPars1
miR1SimPars2[3] <- miR1SimPars2[3]*2
miR1SimPars2[4] <- miR1SimPars2[4]*2

#minus miR with actd
miR1SimPars3 <- miR1SimPars1
miR1SimPars3[2] <- 0

#plus miR with actd
miR1SimPars4 <- miR1SimPars2
miR1SimPars4[2] <- 0

simDatNomiR <- as_tibble(Simulation(miR1SimPars1,1:6000))
simDat2 <- as_tibble(Simulation(miR1SimPars2,1:6000))

initial_state = as.numeric(simDatNomiR[5,])[-1]


simDatNomiRActD <- as_tibble(SimulationActD(miR1SimPars3, 1:6000, initial_state*10))
simDat4 <- as_tibble(SimulationActD(miR1SimPars4, 1:6000, initial_state*10))

simMeans <- tibble(time = simDatNomiR$time,
                    nomiR = apply(simDatNomiR[,-1],1,function(x){weighted.mean(250:0,x)}),
                    withmiR = apply(simDat2[,-1],1,function(x){weighted.mean(250:0,x)}),
                    nomiRSum = rowSums(simDatNomiR[,-1]),
                    withmiRSum = rowSums(simDat2[,-1]))

simMeans$ratio <- log2(simMeans$withmiR/simMeans$nomiR)
simMeans$ratioSum <- log2(simMeans$withmiRSum/simMeans$nomiRSum)


simMeansGather <- gather(simMeans,key = "paramSet", value = "simMean",-time)
simMeansGather$type = "sim"
simMeansGather[which(simMeansGather$paramSet %in% c("nomiRSum","withmiRSum")),]$type = "sum"
simMeansGather[which(simMeansGather$paramSet == "ratio"),]$type = "ratio"
simMeansGather[which(simMeansGather$paramSet == "ratioSum"),]$type = "ratioSum"

simMeansActD <- tibble(time = simDatNomiRActD$time,
                    nomiR = apply(simDatNomiRActD[,-1],1,function(x){weighted.mean(250:0,x)}),
                    withmiR = apply(simDat4[,-1],1,function(x){weighted.mean(250:0,x)}),
                    nomiRSum = rowSums(simDatNomiRActD[,-1]),
                    withmiRSum = rowSums(simDat4[,-1]))

simMeansActD$ratio <- log2(simMeansActD$withmiR/simMeansActD$nomiR)
simMeansActD$ratioSum <- log2(simMeansActD$withmiRSum/simMeansActD$nomiRSum)


simMeansActDGather <- gather(simMeansActD,key = "paramSet", value = "simMean",-time)
simMeansActDGather$type = "sim"
simMeansActDGather[which(simMeansActDGather$paramSet %in% c("nomiRSum","withmiRSum")),]$type = "sum"
simMeansActDGather[which(simMeansActDGather$paramSet == "ratio"),]$type = "ratio"
simMeansActDGather[which(simMeansActDGather$paramSet == "ratioSum"),]$type = "ratioSum"

p1 <- ggplot(simMeansGather[which(simMeansGather$type == "sim"),],aes(x = time, y = simMean,color = paramSet)) + 
    geom_line() + 
    scale_y_continuous(expand = c(0,0),limits =c(0,200),name = "Tail length (nt)") +
    scale_x_continuous(expand = c(0,0),"Time") +
    theme_tim()

OneMRNAtibble <- tibble(time = simDatNomiR$time,
                    nomiR = miR1SimPars1[1] - miR1SimPars1[3]*simDatNomiR$time,
                    withmiR = miR1SimPars2[1] - miR1SimPars2[3]*simDatNomiR$time)

OneMRNAtibbleGat <- gather(OneMRNAtibble,key = "paramSet", value = "simMean",-time)

p2 <- ggplot(simMeansActDGather[which(simMeansActDGather$type == "sim"),],aes(x = time, y = simMean,color = paramSet)) + 
    geom_line() + 
    scale_y_continuous(expand = c(0,0),limits =c(0,200),name = "Tail length (nt)") +
    scale_x_continuous(expand = c(0,0),"Time") +
    theme_tim()

ggsave(p1, file = "figures/other/otherV9/miRSimGraphicalAbstractLeft.pdf",width = 1.75, height = 3)
ggsave(p2, file = "figures/other/otherV9/miRSimGraphicalAbstractRight.pdf",width = 1.75, height = 3)