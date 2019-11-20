
library(deSolve)
library(numDeriv)
library(tictoc)
library(rootSolve)

dyn.load("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/KinetTlnAnalysis/Model/UnLinkV2/AnalyticalTailsPulseStUnlinked.so") #load the model dynamically

Simulation<-function(pars,time){ ##cleanup
  N = 1
  st       = pars[1]
  a        = pars[2]
  k        = pars[3]
  b        = pars[4]
  size     = 12.66052 #These parameters are from fitting with the steady state tail seq data. 
  location = 270.37422
  scale    = 10.80678
  timeSS = 6000
  parameters = c(st,a,k,b,size,location,scale)
  max_tail = 251
  initial_state <- rep(0,max_tail*N)
  tails <- ode.band(func = "ode_deriv", y = initial_state, parms = parameters, times=c(0,time),method = "lsode",bandup = 0, banddown = 1,nspec = max_tail*N, dllname = "AnalyticalTailsPulseStUnlinked",nout=1,initfunc = "ode_p_init",hmax = 60)
  tailsSS <- steady(func = "ode_deriv", y = runif(max_tail*N), parms = parameters, times=c(timeSS),method = "stode",bandup = 0, banddown = 1,nspec = max_tail*N, dllname = "AnalyticalTailsPulseStUnlinked",nout=1,initfunc = "ode_p_init")$y
  tails <- rbind(tails[,-c(1,252)],tailsSS)
  return(tails)

}
# break
# fracSim <- function(x,y,z){
#   pars <- c(x,1E-4,y,z)
#   time <- 1:6000
#   sim <- Simulation(pars,time)
#   frac <- sim[6000,250]/sum(sim[6000,])
#   return(frac)}
# x = seq(50,200,10)
# z = 10^(seq(-3,3,0.1))
# fracDF = data.frame(x = x)
# fracDF$sim <- apply(fracDF,1,fracSim)

SimulationSS<-function(pars,time){ ##cleanup
  N = 1
  st       = pars[1]
  a        = pars[2]
  k        = pars[3]
  b        = pars[4]
  size     = 12.66052 #These parameters are from fitting with the steady state tail seq data. 
  location = 270.37422
  scale    = 10.80678
  time = time[length(time)]
  print(time)
  parameters = c(st,a,k,b,size,location,scale)
  max_tail = 251
  initial_state <- runif(max_tail*N)
  tailsSS <- steady(func = "ode_deriv", y = initial_state, parms = parameters, times=c(time),method = "stode",bandup = 0, banddown = 1,nspec = max_tail*N, dllname = "AnalyticalTailsPulseStUnlinked",nout=1,initfunc = "ode_p_init")
  return(tailsSS$y)

}
