#Figure 4B, the model fitting figure. Attempting to institute a lower baseline.  
#Some notes on this figure from 2019 03 08. The ceiling, floor, and round functions make the script pretty sensitive to the quantile values.
#I decided to just use the residuals, not the Rsq values, for the quantiles. 


library(tictoc)
library(deSolve)
library(numDeriv)
library(reshape2)
library(plyr)
library(DescTools)
library(tibble)
library(tidyverse)
library(ggplus)
# system("R CMD SHLIB /lab/solexa_bartel/teisen/RNAseq/Scripts/models/compiled/GlobalODE_6_plogis_linked.c") #DON'T COMPILE THIS HERE
dyn.load("/lab/solexa_bartel/teisen/RNAseq/Scripts/models/compiled/GlobalODE_6_plogis_linked.so") #load the model dynamically
dyn.load("/lab/solexa_bartel/teisen/RNAseq/Scripts/models/compiled/GlobalODE_6_plogis.so") #load the model dynamically
dyn.load("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/KinetTlnAnalysis/Model/UnLink/AnalyticalTailsPulseStUnlinked.so")
options(warn = 1)
args<-commandArgs(trailingOnly=TRUE)

source("/lab/solexa_bartel/teisen/RNAseq/Scripts/general/ggplot_theme.R")

annot<-read.table("/lab/solexa_bartel/teisen/RNAseq/Annotation_files/TAIL_seq/symbol_to_accession.txt",head=TRUE)

# Simulation<-function(pars,time){ ##cleanup

#   st       = pars[1]
#   a        = pars[2]
#   k        = pars[3]
#   b        = 0.965050588
#   size     = 12.80082479
#   location = 254.7199822
#   scale    = 10.38658404

#   parameters = c(st,a,k,b,size,location,scale)

#   N = 1
#   max_tail = 251
#   initial_state <- rep(0,max_tail*N)
#   tails <- ode.band(func = "ode_deriv", y = initial_state, parms = parameters, times=c(0,time),method = "lsode",bandup = 0, banddown = 1,nspec = max_tail*N, dllname = "GlobalODE_6_plogis_linked",nout=1,initfunc = "ode_p_init",hmax = 60)

#   columns_to_remove = c(1,2,ncol(tails))
#   sim <- tails[-1,-columns_to_remove]
#   return(data.frame(t(sim)))

# }

SimulationUL<-function(pars,time){ ##cleanup

  st       = pars[1]
  a        = pars[2]
  k        = pars[3]
  b        = pars[4]
  size     = 16.16483 #17.27453 #Global fitted parameters as of 2019 07 15. 20 datasets total.
  location = 262.96905 #264.97334
  scale    = 11.07612 #11.08839 

  parameters = c(st,a,k,b,size,location,scale)

  N = 1
  max_tail = 251
  initial_state <- rep(0,max_tail*N)
  tails <- ode.band(func = "ode_deriv", y = initial_state, parms = parameters, times=c(0,time),method = "lsode",bandup = 0, banddown = 1,nspec = max_tail*N, dllname = "AnalyticalTailsPulseStUnlinked",nout=1,initfunc = "ode_p_init",hmax = 0.1,maxsteps=5000000)
  print(dim(tails))
  columns_to_remove = c(1,2,(max_tail - 18):ncol(tails))
  #Note that these values are all replaced by the mean here for the last 20 nt.
  last20ntMean = rep(mean(tails[7,(max_tail - 18):(ncol(tails)-1)]), 20)
  sim <- tails[-1,-columns_to_remove]
  #Add the last 20 nt back to the flattened array. 
  sim <- c(c(t(sim)),last20ntMean)
  return(sim)

}


##Uncomment below to run 4 gene-specific constant dataset
pulse_rate_constants<-read.table("rate_constant_measurements/miR-155_minus_samples_UNLINKV3_V81H5_run3_high_precision_reformat.txt")

colnames(pulse_rate_constants)<-c("accession","st","a","k","b","residual")
pulse_rate_constants<-add_column(pulse_rate_constants,size = rep(16.25643,nrow(pulse_rate_constants)),.after = "b")
pulse_rate_constants<-add_column(pulse_rate_constants,location = rep(263.95156,nrow(pulse_rate_constants)),.after = "size")
pulse_rate_constants<-add_column(pulse_rate_constants,scale = rep(11.05133,nrow(pulse_rate_constants)),.after = "location")
ULPAR = "UL"
##

colnames(pulse_rate_constants)<-c("accession","st_m","a_m","k_m","b_m","size","location","scale","r_m")
data<-read.table("model_input_files/miR-155_minus_sample_mean_tails_50tags_with_expr_background_subtracted_v7_st_last20rm_HYBRID20190731.txt",head=TRUE)
vars <- c(      #miR-155 minus
  7.426524e-16,
  1.287242e-15,
  3.640165e-15,
  9.202623e-15,
  4.320013e-14,
  2.808383e-13,
  2.808383e-13/6) #The variance of the last value is reduced 6 fold
## This increases the weighting value of these points 6 fold
#Updated variances as of 2019 07 15.


# vars <- rep(1,6)

rss <- data.frame(
  "accession" = data[,1],
  "rss" = apply(data[,-1],1,function(x){
    sqs<-(x-mean(x))^2
    sqs[1:230]     <- sqs[1:230]    / vars[1]
    sqs[231:460]   <- sqs[231:460]  / vars[2]
    sqs[461:690]   <- sqs[461:690]  / vars[3]
    sqs[691:920]   <- sqs[691:920]  / vars[4]
    sqs[921:1150]  <- sqs[921:1150] / vars[5]
    sqs[1151:1380] <- sqs[1151:1380]/ vars[6]
    sqs[1381:1400] <- sqs[1381:1400]/ vars[7] #Added in V4 by TJE on 2019 04 01. 
    rss = sum(sqs*1E6)
    })
  )

print(head(rss))

m<-merge(pulse_rate_constants,data,by="accession")

m<-m[complete.cases(m),]
m<-merge(m,rss,by="accession")
m$r_m<-1-m$r_m/m$rss
# m$r_m#<-m$rss/1E6
m$rss<-NULL
# m<-m[which(m$r_m>0.9),] #cutoff?
printNo = 4 
r_l<-NULL
print(quantile(m$r_m))
print(nrow(m))
m<-m[order(m$r_m),]
print(head(m$r_m))

# m <- m[which(!m$accession=="NM_001256521"),]
# quantile_genes<-round(nrow(m)*c(341/1000,500/1000,841/1000))
# quantile_genes<-round(nrow(m)*c(1/4,1/2,3/4))
# quantile_genes<-ceiling(nrow(m)*c(1/5,2/5,3/5,4/5))
# quantile_genes<-floor(nrow(m)*c(1/10,1/4,1/5,3/4,9/10))
# quantile_genes<-round(nrow(m)*c(1:(printNo-1)/printNo))

#Floor or ceiling? Or Round?
quantile_genes<-floor(nrow(m)*c(1/10,1/4,3/4,9/10))


# quantile_genes<-which(m$accession %in% c("NM_007907","NM_007393","NM_009716"))
# quantile_genes <- which(m$accession %in% c("NM_010833","NM_172145","NM_134131","NM_026267"))
# quantile_genes <- which(m$accession %in% c("NM_010849","NM_008037","NM_010499","NM_010500")) #IEGs
print(quantile_genes)

# quantile_genes <- which(m$accession %in% c("NM_011897",
# "NM_001033457",
# "NM_010143",
# "NM_025825",
# "NM_178633",
# "NM_178798",
# "NM_009302",
# "NM_011519",
# "NM_019773",
# "NM_025380",
# "NM_026975",
# "NM_026179",
# "NM_025550",
# "NM_008252",
# "NM_053104",
# "NM_181392",
# "NM_026441",
# "NM_008212",
# "NM_027959")
# )

p1 <- ggplot(m,aes(x=r_m)) +
geom_histogram(bins=100) + 
theme_tim() +
scale_y_continuous(name=NULL,expand=c(0,0)) +
scale_x_continuous(expand=c(0,0)) +
geom_vline(linetype="dashed",color="black",xintercept = m[quantile_genes,]$r_m)
# write.table(m,file="RsquaredValues.txt",row.names=FALSE,sep="\t",quote=FALSE)
ggsave(plot=p1,file="figures/version_9/FigS7A.pdf",width=2,height=2,useDingbats=FALSE)

bins <- 2 #How many nucleotides to bin?
bin_vector<-function(bins,vec){
    val <- 1:(length(vec)/bins)
    unlist(lapply(val,function(x,vec){rep(sum(vec[(bins*(x-1)+1):(bins*x)]),bins)},vec=vec))} #Sum

tmeans = NULL
for(gene in quantile_genes){
  pulse_params<-as.numeric(m[gene,2:8])
  if(ULPAR=="UL"){sim<-SimulationUL(pulse_params,time=c(40,60,120,240,480,6000)-35)}

  t1 <- data.frame(
    c(sim[1:230]   ,rep(0,20)),
    c(sim[231:460] ,rep(0,20)),
    c(sim[461:690] ,rep(0,20)),
    c(sim[691:920] ,rep(0,20)),
    c(sim[921:1150],rep(0,20)),
    sim[1151:1400]
    )

  print(dim(t1))
  colnames(t1)<-c(40,60,120,240,480,6000)-35 
  ssMean20 =  mean(as.numeric(m[gene,1390:1409]))
  rd <- data.frame(
    c(t(m[gene,10:239]),rep(0,20)),
    c(t(m[gene,240:469]),rep(0,20)),
    c(t(m[gene,470:699]),rep(0,20)),
    c(t(m[gene,700:929]),rep(0,20)),
    c(t(m[gene,930:1159]),rep(0,20)),
    c(t(m[gene,1160:1389]),rep(ssMean20,20)))

  colnames(rd)<-c(40,60,120,240,480,6000)-35  

  t1$sample<-"model"
  rd$sample<-"data"



  t1$tl<-249:0
  rd$tl<-249:0

  t1m <- data.frame(apply(t1[,1:6],2,function(x){weighted.mean(t1$tl,w = x)}))
  r1m <- data.frame(apply(rd[,1:6],2,function(x){weighted.mean(rd$tl,w = x)}))


  tmeansTemp = data.frame(
    tp = c(40,60,120,240,480,6000)-35,
    model = t1m,
    data = r1m
    )

  tmeansTemp$gene = m[gene,1]

  tmeans <- rbind(tmeans,tmeansTemp)


  t2l<-rbind(t1,rd)
  t2l$gene<-m[gene,1]
  print(m[gene,1])
  print(c("Residual",m[gene,9]),sep="\t",quote=FALSE)
  print("-----------------------------")
  t2l$residual <- m[gene,9]
  r_l<-rbind(r_l,t2l)
}


r_l$gene<-factor(r_l$gene,levels=unique(r_l[order(r_l$residual,decreasing=FALSE),]$gene))

r_l_melt<-melt(r_l,id=c("gene","tl","sample","residual"))
# r_l_melt$smooth<-unlist(lapply(1:(nrow(r_l_melt)/bins),bin_vector))
cols=c("black","red","blue")

colnames(tmeans) = c("time","model","data","accession")
# r_l_melt$ymax = 1E-6

# r_l_melt[which(r_l_melt$gene=="NM_139063"),]$ymax = 1.5E-7
# r_l_melt[which(r_l_melt$gene=="NM_001199275"),]$ymax = 3E-7
# r_l_melt[which(r_l_melt$gene=="NM_009396"),]$ymax = 5E-7
# r_l_melt[which(r_l_melt$gene=="NM_026932"),]$ymax = 7.5E-7

colnames(annot)[2]<- "gene"
r_l_melt <- merge(r_l_melt,annot,by="gene")

r_l_melt$symbol<-factor(r_l_melt$symbol,levels=unique(r_l_melt[order(r_l_melt$residual,decreasing=FALSE),]$symbol))

print(head(r_l_melt))
r_l_melt$value <- r_l_melt$value * 1E7
r_l_melt <- r_l_melt[order(r_l_melt$tl),]
r_l_melt <- r_l_melt[order(r_l_melt$variable),]
r_l_melt <- r_l_melt[order(r_l_melt$symbol),]
r_l_melt <- r_l_melt[order(r_l_melt$sample),]


r_l_melt$names <- r_l_melt$symbol
levels(r_l_melt$variable) <- c("40 min","1 h","2 h","4 h","8 h","Ss")
levels(r_l_melt$symbol) <- rev(c("10th percentile","25th percentile","75th percentile","90th percentile"))

# r_l_melt <- r_l_melt[which(!r_l_melt$tl %in% 0:7),] #comment this line and change the one below to 0 to display the LIM1. 
LIM1 = 0
msg <- r_l_melt[-which(duplicated(r_l_melt$gene)),]


msg$yend = c(
  max(bin_vector(bins,r_l_melt[r_l_melt$gene == msg[1,]$gene,]$value)),
  max(bin_vector(bins,r_l_melt[r_l_melt$gene == msg[2,]$gene,]$value)),
  max(bin_vector(bins,r_l_melt[r_l_melt$gene == msg[3,]$gene,]$value)),
  max(bin_vector(bins,r_l_melt[r_l_melt$gene == msg[4,]$gene,]$value))
  # max(r_l_melt[r_l_melt$gene == msg[5,]$gene,]$value)

  )

msgX <- distinct( r_l_melt,gene,variable,.keep_all = TRUE)


msgX$xbegin = rep(c(
  20,
  20,
  20,
  20,
  20,
  0
  ),4)

print(head(msg))

r_l_melt$symbol <- factor(r_l_melt$symbol,levels = unique(r_l_melt[order(-r_l_melt$residual),]$symbol))

r_l_melt_mod <- r_l_melt[which(!(r_l_melt$tl < 20 & r_l_melt$variable %in% c("40 min","1 h","2 h","4 h","8 h"))),]
# r_l_melt_mod[r_l_melt_mod<0] <- 0
p8<-ggplot()+
  scale_y_continuous(name=NULL,expand=c(0,0))+
  scale_x_continuous(name=NULL,expand=c(0,0),breaks = c(0,20,50,100,150,200,250),labels = c(0,20,50,100,150,200,250))+
  # geom_step(data = r_l_melt[which(r_l_melt$sample == "data"),],aes(x=tl,y=value),size=0.2,color="red")+

  geom_step(data = r_l_melt_mod[which(r_l_melt_mod$sample == "data"),],aes(x=tl,y=bin_vector(bins,value)),size=0.2,color="red")+
  # geom_smooth(data = r_l_melt_mod[which(r_l_melt_mod$sample == "data"),],aes(x=tl,y=bin_vector(bins,value)),size=0.5,method="loess",color="black",span=0.75)+ #Smoothing kernel
  # geom_line(data = r_l_melt_mod[which(r_l_melt_mod$sample == "model"),],aes(x=tl,y=bin_vector(bins, value)),size=0.6,color="blue",alpha = 0.8)+
  geom_line(data = r_l_melt_mod[which(r_l_melt_mod$sample == "model"),],aes(x=tl,y = value * bins), size=0.6,color="blue",alpha = 0.8)+
  theme_tim()+
  facet_grid(symbol~variable,scale="free")+
  geom_segment(data=msg,aes(x=0,xend=0,y=0,yend=yend),size=0.5)+
  geom_segment(data=msgX,aes(x=xbegin,xend=250,y=0,yend=0),size=0.5)+
  geom_text(data=msg,aes(label=names,x=80,y=yend/2),size=6*5/14)+
  theme(axis.line.x=element_blank(),
        axis.line.y=element_blank())

# if(ULPAR=="ULGL"){ggsave(plot=p8,file="/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/figures/other/Fig4CtempUnlinkedGlobal.pdf",width=6.5,height=printNo,useDingbats=FALSE)}
if(ULPAR=="UL"){ggsave(plot=p8,file="/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/figures/version_9/Fig4C.pdf",width=5.5,height=printNo,useDingbats=FALSE)}
# if(ULPAR=="GL"){ggsave(plot=p8,file="/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/figures/other/Fig4CtempLinked.pdf",width=6.5,height=printNo,useDingbats=FALSE)}
# ggsave(plot=p2,file="figures/version_5/fig5D_panel_2_v2.pdf",width=3.5,height=1.5,useDingbats=FALSE)
