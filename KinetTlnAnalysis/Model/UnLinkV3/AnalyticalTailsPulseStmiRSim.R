#populate vector of tails

library(deSolve)
library(numDeriv)
library(tidyverse)
system("R CMD SHLIB /lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/KinetTlnAnalysis/Model/UnLinkV3/AnalyticalTailsPulseStmiRSim.c")
dyn.load("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/KinetTlnAnalysis/Model/UnLinkV3/AnalyticalTailsPulseStmiRSim.so") #load the model dynamically
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
  locationmir    = pars[8] 
  scalemir       = pars[9] 
  bmir           = pars[10] 

  parameters = c(st,a,k,b,size,location,scale,locationmir,scalemir,bmir)
  max_tail = 251

  initial_state <- rep(0,max_tail) #All abundances begin with 0. 
  #The simulation, passed to c code called ode_deriv, using lsode.
  #This is for a banded jacobian.
  #Note hmax has a major impact on memory usage, time, and precision. 
  tails <- ode.band(func = "ode_deriv", y = initial_state, parms = parameters, 
          times = c(0,time), method = "lsode",bandup = 0, banddown = 1,
          nspec = max_tail, dllname = "AnalyticalTailsPulseStmiRSim",
          nout=1, initfunc = "ode_p_init", hmax = 1,maxsteps=5000000)

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
  0.05,
  0.1,
  50,
  260,
  10,
  180,
  10,
  0.0
  )

miR1SimPars2 <- miR1SimPars1
miR1SimPars2[10] <- 0.005


# miR1SimPars2a <- miR1SimPars2
# miR1SimPars2b <- miR1SimPars2
# miR1SimPars2c <- miR1SimPars2
# miR1SimPars2d <- miR1SimPars2
# miR1SimPars2e <- miR1SimPars2

# miR1SimPars2a[3] <- 0.001
# miR1SimPars2b[3] <- 0.01
# miR1SimPars2c[3] <- 0.1
# miR1SimPars2d[3] <- 1
# miR1SimPars2e[3] <- 10


# miRSimParsAll <- list(
#   miR1SimPars2a,
#   miR1SimPars2b,
#   miR1SimPars2c,
#   miR1SimPars2d,
#   miR1SimPars2e)




simDatNomiR <- as_tibble(Simulation(miR1SimPars1,1:6000))
simDat2 <- as_tibble(Simulation(miR1SimPars2,1:6000))

simDatNomiRSs <- simDatNomiR[6000,]
simDatmiRSs <- simDat2[6000,]

# simDatNomiRSs[,-1] <- simDatNomiRSs[,-1]/sum(simDatNomiRSs[,-1])
# simDatmiRSs[,-1] <- simDatmiRSs[,-1]/sum(simDatmiRSs[,-1])

allCDF <- tibble(
  length = rep(0:250,2),
  mir = rep(c("no_mir","mir"),each = 251),
  cdfVec = c(cdff(simDatNomiRSs[,-1]),cdff(simDatmiRSs[,-1]))
  )

pCdf <- ggplot(allCDF,aes(x = length, y = cdfVec, color = mir)) +
  geom_line() +
  scale_x_continuous(expand = c(0,0),name = "Tail length (nt)" )+
  scale_y_continuous(expand = c(0,0),name = "Cumulative abundance") +
  theme_tim()

ggsave(pCdf, file = "figures/other/otherV8/miR1_longTailSim_CDF.pdf",width = 2.5, height = 2.5)


# simGather <- gather(simDatNomiR,key = "tail_length",value = "abundance",-time)
# simGather2 <- gather(simDat2,key = "tail_length",value = "abundance",-time)
# simGather$paramSet <- "nomiR"
# simGather2$paramSet <- "withmiR"

# simall <- bind_rows(simGather,simGather2)
# simall$tail_length <- as.numeric(simall$tail_length)

# simDat2a <- as_tibble(Simulation(miR1SimPars2a,1:6000))
# simDat2b <- as_tibble(Simulation(miR1SimPars2b,1:6000))
# simDat2c <- as_tibble(Simulation(miR1SimPars2c,1:6000))
# simDat2d <- as_tibble(Simulation(miR1SimPars2d,1:6000))
# simDat2e <- as_tibble(Simulation(miR1SimPars2e,1:6000))

# simDat2a$data <- 0.001
# simDat2b$data <- 0.01
# simDat2c$data <- 0.1
# simDat2d$data <- 1
# simDat2e$data <- 10

# simallVars <- bind_rows(
#   simDat2a,
#   simDat2b,
#   simDat2c,
#   simDat2d,
#   simDat2e)

# simallVars$mean <- apply(select(simallVars,-time,-data),1,function(x){weighted.mean(250:0,x)})
# simallVars$sum <- apply(select(simallVars,-time,-data,-mean),1,sum)

# simallVars <- select(simallVars,time,data,mean,sum)

# simallVarsMerge <- mutate(simDatNomiR,meanNomiR = apply(select(simDatNomiR,-time),1,function(x){weighted.mean(250:0,x)})) %>% mutate(rs = rowSums(select(simDatNomiR,-time))) %>% select(time,meanNomiR,rs) %>% right_join(simallVars,by = "time")


# # txFunc <- tibble(
# #         length = 250:0,
# #         nomiR = dnbinom(250:0,mu = 180,50)*1E-7,
# #         withmiR = dnbinom(250:0,mu = 180,50)*1E-7)

# # txFuncGather <- gather(txFunc,key = "paramSet", value = "txFunc", -length)

# # dcpFunc <- tibble(
# #         length = 250:0,
# #         nomiR = plogis(0:250,260,10)*0.1,
# #         withmiR = plogis(0:250,260,10)*0.1 + plogis(250:0,160,10)*0.005)

# # dcpFuncGather <- gather(dcpFunc,key = "paramSet", value = "dcpFunc", -length)

# simMeans <- tibble(time = simDatNomiR$time,
#                     nomiR = apply(simDatNomiR[,-1],1,function(x){weighted.mean(250:0,x)}),
#                     withmiR = apply(simDat2[,-1],1,function(x){weighted.mean(250:0,x)}),
#                     nomiRSum = rowSums(simDatNomiR[,-1]),
#                     withmiRSum = rowSums(simDat2[,-1]))

# simMeans$ratio <- log2(simMeans$withmiR/simMeans$nomiR)
# simMeans$ratioSum <- log2(simMeans$withmiRSum/simMeans$nomiRSum)


# simMeansGather <- gather(simMeans,key = "paramSet", value = "simMean",-time)
# simMeansGather$type = "sim"
# simMeansGather[which(simMeansGather$paramSet %in% c("nomiRSum","withmiRSum")),]$type = "sum"
# simMeansGather[which(simMeansGather$paramSet == "ratio"),]$type = "ratio"
# simMeansGather[which(simMeansGather$paramSet == "ratioSum"),]$type = "ratioSum"


# # p1 <- ggplot(simMeansGather[which(simMeansGather$type == "sim"),],aes(x = time, y = simMean,color = paramSet)) + 
# #     geom_line() + 
# #     scale_y_continuous(expand = c(0,0),limits =c(50,200),name = "Tail length (nt)") +
# #     scale_x_continuous(expand = c(0,0),"Time") +
# #     theme_tim_label()

# # p1a <- ggplot(simMeansGather[which(simMeansGather$type == "ratio"),],aes(x = time, y = simMean)) + 
# #     geom_line() + 
# #     scale_y_continuous(expand = c(0,0),name = "Mean fold change (log2)") +
# #     scale_x_continuous(expand = c(0,0),name = "Time") +
# #     theme_tim_label()

# # p1b <- ggplot(simMeansGather[which(simMeansGather$type == "ratioSum"),],aes(x = time, y = simMean)) + 
# #     geom_line() + 
# #     scale_y_continuous(expand = c(0,0),name = "Mean fold change (log2)") +
# #     scale_x_continuous(expand = c(0,0),name = "Time") +
# #     theme_tim_label()

# # p2 <- ggplot(simallVars,aes(x = time, y = mean,,color = as.factor(data))) + 
# #     geom_line() + 
# #     scale_y_continuous(expand = c(0,0),limits =c(0,200),name = "Tail length (nt)") +
# #     scale_x_continuous(expand = c(0,0),"Time") +
# #     theme_tim_label()

# # p2a <- ggplot(simallVarsMerge,aes(x = time, y = log2(mean/meanNomiR),color = as.factor(data))) + 
# #     geom_line() + 
# #     scale_y_continuous(expand = c(0,0),name = "Mean fold change (log2)") +
# #     scale_x_continuous(expand = c(0,0),name = "Time") +
# #     theme_tim_label()

# # p2b <- ggplot(simallVarsMerge,aes(x = time, y = log2(sum/rs),,color = as.factor(data))) + 
# #     geom_line() + 
# #     scale_y_continuous(expand = c(0,0),name = "Mean fold change (log2)") +
# #     scale_x_continuous(expand = c(0,0),name = "Time") +
# #     theme_tim_label()

# # p2 <- ggplot(dcpFuncGather,aes(x = length, y = dcpFunc,color = paramSet)) + 
# #     geom_line() + 
# #     scale_y_continuous(expand = c(0,0)) +
# #     scale_x_continuous(expand = c(0,0)) +
# #     theme_tim()

# # p3 <- ggplot(txFuncGather,aes(x = length, y = txFunc,color = paramSet)) + 
# #     geom_line() + 
# #     scale_y_continuous(expand = c(0,0)) +
# #     scale_x_continuous(expand = c(0,0)) +
# #     theme_tim()

# # p4 <- ggplot(simall[which(simall$time == 5500),],aes(x = tail_length, y = abundance, color = paramSet)) + 
# #   geom_line() +
# #   scale_y_continuous(expand = c(0,0)) +
# #   scale_x_continuous(expand = c(0,0)) +

# #   theme_tim()



# # ggsave(p1, file = "figures/other/otherV8/miR1_longTailSim_means.pdf",width = 2.5, height = 2.5)
# # ggsave(p1a, file = "figures/other/otherV8/miR1_longTailSim_ratio.pdf",width = 2.5, height = 2.5)
# # ggsave(p1b, file = "figures/other/otherV8/miR1_longTailSim_sum.pdf",width = 2.5, height = 2.5)

# # ggsave(p2, file = "figures/other/otherV8/miR1_longTailSim_means_multi.pdf",width = 2.5, height = 2.5)
# # ggsave(p2a, file = "figures/other/otherV8/miR1_longTailSim_ratio_multi.pdf",width = 2.5, height = 2.5)
# # ggsave(p2b, file = "figures/other/otherV8/miR1_longTailSim_sum_multi.pdf",width = 2.5, height = 2.5)

# # ggsave(p2, file = "figures/other/otherV8/miR1_longTailSim_dcp_func.pdf",width = 3, height = 3)
# # ggsave(p3, file = "figures/other/otherV8/miR1_longTailSim_tx_func.pdf",width = 3, height = 3)
# # ggsave(p4, file = "figures/other/otherV8/miR1_longTailSim_Ss_distr.pdf",width = 3, height = 3)