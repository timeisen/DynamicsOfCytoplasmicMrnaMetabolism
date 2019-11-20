#Fig1G
library(ggplot2)
library(data.table)
options("scipen"=10)

#Analysis of Chang et al., 2014

mirdata<-read.table("other_analyses/Changetal2014Data/HeLaFcData20190819_data_compiled.txt",head=TRUE,sep="\t")
mirdata$Data<-factor(mirdata$Data,levels=c("RNA","PAL"))
source("ggplot_theme.R")

fill_cols_miR1<-c("#F8766D","#00BA38")




p1<-ggplot(mirdata[which(mirdata$microRNA=="miR1"),],aes(x=Data,y=tpUTR,group=Time.Point,fill=Data))+
geom_bar(width=0.6,position = position_dodge(width=0.8),stat='identity')+
scale_y_continuous(limits=c(-1,0.2),breaks=round(seq(-0.3,.0,length.out=6),1),name=NULL,expand = c(0,0))+
geom_segment(aes(y=-0.3,yend=0.0,x=-Inf,xend=-Inf))+theme_tim()+theme(axis.line = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())+
geom_errorbar(data=mirdata[which(mirdata$microRNA=="miR1"),],aes(ymin=tpUTR-tpUTR_err, ymax=tpUTR+tpUTR_err),
	width=.2,                    # Width of the error bars
	position=position_dodge(.8))+geom_text(data=mirdata[which(mirdata$microRNA=="miR1"),],aes(label=Site.asteriks,y=tpUTR-0.1),position = position_dodge(width=0.8),vjust=0) +
scale_fill_manual(values = fill_cols_miR1)

p2<-ggplot(mirdata[which(mirdata$microRNA=="miR1"),],aes(x=Data,y=X25pct,group=Time.Point,fill=Data))+
geom_bar(width=0.6,position = position_dodge(width=0.8),stat='identity')+
scale_y_continuous(limits=c(-1,0.2),breaks=seq(-.6,.0,by=0.1),name=NULL,expand = c(0,0))+
geom_segment(aes(y=-0.6,yend=0.0,x=-Inf,xend=-Inf))+
theme_tim()+theme(axis.line = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())+
geom_errorbar(data=mirdata[which(mirdata$microRNA=="miR1"),],aes(ymin=X25pct-X25pct_err, ymax=X25pct+X25pct_err),
	width=.2,                    # Width of the error bars
	position=position_dodge(.8))+geom_text(data=mirdata[which(mirdata$microRNA=="miR1"),],aes(label=Top.site.asteriks,y=X25pct-0.15),position = position_dodge(width=0.8),vjust=0) +
scale_fill_manual(values = fill_cols_miR1)

ggsave(plot=p1,file="/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/figures/version_9/fig1EmiRNAleft.pdf",width = 1.5,height = 2.33)
ggsave(plot=p2,file="/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/figures/version_9/fig1EmiRNAright.pdf",width = 1.5,height = 2.33)

