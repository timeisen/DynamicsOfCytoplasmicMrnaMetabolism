#Produces a plot of mean tail lengths over time.
library(ggplot2)
library(reshape2)
source("ggplot_theme.R")

m <- read.table("processed_files/miR-155_minus_sample_mean_tails_50tags_background_subtracted_v7.txt")
colnames(m) <- c("accession","tp40m","tp1h","tp2h","tp4h","tp8h","tpSS")

mQ <- data.frame(t(apply(m[,-1],2,function(x){quantile(x,probs=c(0,0.05,0.25,0.5,0.75,0.95))})))
mQ$time <- c(40,60,120,240,480,1000)

mm <- melt(m,id="accession")
mm$time = NA
mm[which(mm$variable=="tp40m"),]$time = 40
mm[which(mm$variable=="tp1h"),]$time = 60
mm[which(mm$variable=="tp2h"),]$time = 120
mm[which(mm$variable=="tp4h"),]$time = 240
mm[which(mm$variable=="tp8h"),]$time = 480
mm[which(mm$variable=="tpSS"),]$time = 1000


## The next two lines remove the steady-state lines. 2019 05 20. 
# mm[which(mm$variable=="tpSS"),]$value = median(mm[which(mm$variable=="tpSS"),]$value)
# mm[which(mm$variable=="tpSS"),]$accession = NA


mm[which(mm$value > 250),]$value = 250
mm[which(mm$value < 50),]$value = 50

#Ignore the bizarre warning about width. 
p1 <- ggplot()+
geom_line(data = mm,aes(x=time,y=value,group=accession),alpha=0.01,size = 0.2)+
geom_boxplot(data=mQ,aes(x=time,lower=X25.,middle=X50.,upper = X75.,ymin=X5.,ymax=X95.,group=time),stat="identity",alpha=0.8,size=0.2, position = position_dodge2(preserve = "total"), width = 0.1,color = "#1b75bc")+
scale_y_continuous(limits=c(50,250),expand=c(0,0))+
scale_x_continuous(trans="log10",breaks=c(40,60,120,240,480,1000),labels=c("0.7","1","2","4","8","Ss"))+
theme_tim()+
theme(axis.line.x = element_blank())+
geom_segment(aes(x = 0,xend = 480,y = 50, yend = 50))

ggsave(plot=p1,file="figures/version_9/FigureSXMeanConnections.pdf",width = 50.8,height = 31.75,useDingbats=FALSE, units = 'mm')