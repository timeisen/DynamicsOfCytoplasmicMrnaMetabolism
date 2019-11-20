#miRNA fig 6
library(ggplot2)
library(data.table)
options("scipen"=10)
source("ggplot_theme.R")

mirdata<-fread("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/miR_analysis/fc_simulation/rateComparisonsSsFitV2/simsV2/RateDataV6.txt",head=TRUE,sep="\t")
mirdata$Data<-factor(mirdata$Data,levels=c("RNA","PAL"))
mirdataraw<-read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/miR_analysis/fig6_updated/miR-data_fig6_updated.txt",head=TRUE,sep="\t")
mirdataraw$Time.Point<-factor(mirdataraw$Time.Point,levels=c("40min","1hr","2hr","4hr","8hr","Steady State"))

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
cols = gg_color_hue(3)


#create a new df with just the time points for comparison
mirdatasim <- mirdata[which(mirdata$Time_point %in% c(40,60,120,240,480,6000)),]
mirdatasim$Time.Point = "Steady State"
mirdatasim[which(mirdatasim$Time_point == 40),]$Time.Point = "40min"
mirdatasim[which(mirdatasim$Time_point == 60),]$Time.Point = "1hr"
mirdatasim[which(mirdatasim$Time_point == 120),]$Time.Point = "2hr"
mirdatasim[which(mirdatasim$Time_point == 240),]$Time.Point = "4hr"
mirdatasim[which(mirdatasim$Time_point == 480),]$Time.Point = "8hr"
mirdatasim <- mirdatasim[which(!(mirdatasim$microRNA=="miR1" & mirdatasim$Time.Point=="8hr")),]
mirdataraw <- mirdataraw[which(mirdataraw$Data != "RPF"),]
mirdataraw[,c(5:8)] <- NULL
mirdatasim[,c(7)] <- NULL
mirdataraw$Rate = "Mea"
m <- rbind(mirdataraw,mirdatasim)
m$Rate <- factor(m$Rate,levels=c("Dcp","Dea","Two","Mea"))
m$Data <- factor(m$Data,levels=c("RNA","PAL"))

m$tp = 30
m[m$Time.Point=="40min",]$tp = 40
m[m$Time.Point=="1hr",]$tp = 60
m[m$Time.Point=="2hr",]$tp = 120
m[m$Time.Point=="4hr",]$tp = 240
m[m$Time.Point=="8hr",]$tp = 480
m[m$Time.Point=="Steady State",]$tp = 535 #or 800?
levels(m$Time.Point) = c("40 m","1 h","2 h","4 h","8 h", "Ss")
mseg <- head(m,1)


p3 <- ggplot(m,aes(x=tp,y=tpUTR,fill=Rate))+
geom_bar(width=0.1,position=position_dodge(width=0.15),stat="identity")+
facet_wrap(~microRNA+Data,ncol=4)+
theme_tim()+
scale_y_continuous(limits=c(-1,1.0),breaks=c(-0.4,-0.2,0,0.2,0.4),name=NULL,expand = c(0,0))+
scale_x_continuous(trans="log10",limits=c(30,1000),expand=c(0,0))+
geom_segment(data=mseg,aes(y=-0.4,yend=0.4,x=30,xend=30))+theme_tim()+theme(axis.line = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())+
geom_segment(aes(y = 0, yend = 0,x = 30, xend = 600),color="black",size=0.1)+
geom_segment(aes(y = 0, yend = 0,x = 630, xend = 1000),color="black",size=0.1)+
geom_errorbar(data=m,aes(ymin=tpUTR-tpUTR_err, ymax=tpUTR+tpUTR_err),size=0.2,width=.1,position=position_dodge(width=0.15))+
scale_color_manual(values = c(cols,"black"))

p4 <- ggplot(m,aes(x=tp,y=X25pct,color=Rate))+
geom_point(size=0.2,shape=16,position=position_dodge(width=0.15))+
facet_wrap(~microRNA+Data,ncol=4)+
theme_tim()+
scale_y_continuous(limits=c(-1,1.0),breaks=c(-0.8,-0.6,-0.4,-0.2,0,0.2,0.4),name=NULL,expand = c(0,0))+
scale_x_continuous(trans="log10",limits=c(30,1000),expand=c(0,0))+
geom_segment(data=mseg,aes(y=-0.8,yend=0.4,x=30,xend=30))+theme_tim()+theme(axis.line = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())+
geom_segment(aes(y = 0, yend = 0,x = 30, xend = 600),color="black",size=0.1)+
geom_segment(aes(y = 0, yend = 0,x = 630, xend = 1000),color="black",size=0.1)+
geom_errorbar(data=m,aes(ymin=X25pct-X25pct_err, ymax=X25pct+X25pct_err),size=0.2,width=.1,position=position_dodge(width=0.15))+
scale_color_manual(values = c(cols,"black"))

# ggsave(plot=p3,file="figures/version_7/FigS7Dv2.pdf",width=5,height=1.5,useDingbats=FALSE)
# ggsave(plot=p4,file="figures/version_7/FigS7Ev2.pdf",width=5,height=1.5,useDingbats=FALSE)
# break
# Create a new dataframe containing just the steady state data, repeated as a horizontal line
mirdataSS <- mirdata[which(mirdata$Time_point==6000),]
mirdataSS <- mirdataSS[rep(seq_len(nrow(mirdataSS)), each=30),]
mirdataSS$Time_point = rep(520:549,24)

# Y axis data frame
geomYaxis <- head(mirdataSS[which(mirdataSS$microRNA=="miR1"),],1)
colnames(geomYaxis)[1:4] <- c("x","xend","y","yend")
geomYaxis$x = 0
geomYaxis$xend = 0
geomYaxis$y = -0.8
geomYaxis$yend = 0.4


colnames(m)[9] <- "Time_point"
m <- m[which(m$Rate == "Mea"),]
# plotting
maxTime = 500
mirdata[mirdata$Time_point>maxTime,]$Time_point=NA
p1 <- ggplot(mirdata,aes(x = Time_point,y = tpUTR,color = Rate))+
geom_ribbon(aes(ymin = tpUTR - tpUTR_err, ymax = tpUTR + tpUTR_err), alpha = 0.1,linetype=0, size = 0.5)+
geom_line(alpha = 0.8, size = 0.5)+
geom_line(alpha = 0.8, size = 0.5,data=mirdataSS,aes(x = Time_point,y = tpUTR, color = Rate))+
geom_ribbon(data=mirdataSS,aes(ymin = tpUTR - tpUTR_err, ymax = tpUTR + tpUTR_err), alpha = 0.1,linetype=0, size = 0.5)+
geom_point(data = m, aes(x = Time_point, y = tpUTR),color="black",size=0.5)+
geom_errorbar(data = m,aes(x = Time_point, ymin = tpUTR - tpUTR_err, ymax = tpUTR + tpUTR_err),color="black",size=0.5,position=position_dodge(width=0.38))+
facet_wrap(~microRNA+Data,nrow=1)+
scale_x_continuous(limits=c(0,550),expand=c(0,0),breaks=c((0:4)*60*2,535),labels=c((0:4)*2,"Ss"))+
scale_y_continuous(limits=c(-1.0,0.8),expand=c(0,0),breaks=c(-0.4,-0.2,0.0,0.2))+
geom_segment(aes(x=0,xend=maxTime,y=0,yend=0),color="black",size=0.5/4)+
geom_segment(aes(x=520,xend=549,y=0,yend=0),color="black",size=0.5/4)+
geom_segment(data=geomYaxis,aes(x=x,xend=xend,y=-0.4,yend=0.2),color="black",size=0.5)+
geom_segment(aes(x=0,xend=8*60,y=-1.0,yend=-1.0),color="black",size=0.5/4)+
theme_tim()+
theme(axis.line.x = element_blank(),
	axis.line.y = element_blank(),
  strip.background = element_blank(),
  strip.text.x = element_blank()
)

p2 <- ggplot(mirdata,aes(x = Time_point,y = X25pct,color = Rate))+
geom_ribbon(aes(ymin = X25pct - X25pct_err, ymax = X25pct + X25pct_err), alpha = 0.1,linetype=0, size = 0.5)+
geom_line(alpha = 0.8, size = 0.5)+
geom_line(alpha = 0.8, size = 0.5,data=mirdataSS,aes(x = Time_point,y = X25pct, color = Rate))+
geom_ribbon(data=mirdataSS,aes(ymin = X25pct - X25pct_err, ymax = X25pct + X25pct_err), alpha = 0.1,linetype=0, size = 0.5)+
geom_point(data = m, aes(x = Time_point, y = X25pct),color="black",size=0.5)+
geom_errorbar(data = m,aes(x = Time_point, ymin = X25pct - X25pct_err, ymax = X25pct + X25pct_err),color="black",size=0.5,position=position_dodge(width=0.38))+
facet_wrap(~microRNA+Data,nrow=1)+
scale_x_continuous(limits=c(0,550),expand=c(0,0),breaks=c((0:4)*60*2,535),labels=c((0:4)*2,"Ss"))+
scale_y_continuous(limits=c(-1,0.8),breaks=c(-0.8,-0.6,-0.4,-0.2,0,0.2,0.4),name=NULL,expand = c(0,0))+
geom_segment(aes(x=0,xend=maxTime,y=0,yend=0),color="black",size=0.5/4)+
geom_segment(aes(x=520,xend=550,y=0,yend=0),color="black",size=0.5/4)+
geom_segment(data=geomYaxis,aes(x=x,xend=xend,y=y,yend=yend),color="black",size=0.5)+
geom_segment(aes(x=0,xend=8*60,y=-1.0,yend=-1.0),color="black",size=0.5/4)+
theme_tim()+
theme(axis.line.x = element_blank(),
	axis.line.y = element_blank(),
  strip.background = element_blank(),
  strip.text.x = element_blank()
)


ggsave(plot=p1,file="figures/version_9/MirSimFigs/FigS7Bv3.png",width=7,height=4)
ggsave(plot=p2,file="figures/version_9/MirSimFigs/FigS7Cv3.png",width=7,height=4)
