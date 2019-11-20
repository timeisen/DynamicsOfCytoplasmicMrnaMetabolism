#3p end tags
library(ggplot2)
library(data.table)
library(scales)
options("scipen"=-3)
source("/lab/solexa_bartel/teisen/RNAseq/Scripts/general/ggplot_theme.R")


log_ticks<-function(start,end){
  tic_num<-log10(end)-log10(start)+1
  labels<-rep(" ",10*tic_num)
  labels[seq(10,(tic_num*10),10)]<-10^(log10(start):log10(end))*10
  breaks<-c(matrix(1:10)%*%t(matrix(10^(log10(start):log10(end)))))
  return(list(labels,breaks))
}


v1_1<-fread("3p_ends/tag_count_files/total_miR-1_palseq_v1_tags_36M.txt")
v2_1<-fread("3p_ends/tag_count_files/total_miR-1_tags.txt")
v2_155<-fread("3p_ends/tag_count_files/total_miR-155_tags.txt")
rnaseq_1<-fread("3p_ends/tag_count_files/total_miR-1_tags_rnaseq.txt")
rnaseq_155<-fread("3p_ends/tag_count_files/total_miR-155_tags_rnaseq.txt")

colnames(v1_1)<-c("tag","count_v1_1")
colnames(v2_1)<-c("tag","count_v2_1")
colnames(v2_155)<-c("tag","count_v2_155")
colnames(rnaseq_1)<-c("tag","count_rnaseq_1")
colnames(rnaseq_155)<-c("tag","count_rnaseq_155")

tech_comp<-merge(v1_1,v2_1,by="tag",all=TRUE)
rep_comp<-merge(v2_1,v2_155,by="tag",all=TRUE)
rnaseq_comp<-merge(rnaseq_1,rnaseq_155,by="tag",all=TRUE)

tech_comp[is.na(tech_comp)]<-0
rep_comp[is.na(rep_comp)]<-0
rnaseq_comp[is.na(rnaseq_comp)]<-0

tech_comp[,c(2,3)]<- (tech_comp[,c(2,3)]+.1)
rep_comp[,c(2,3)]<-rep_comp[,c(2,3)]+.1
rnaseq_comp[,c(2,3)]<-rnaseq_comp[,c(2,3)]+.1


# tech_comp_cutoff<-tech_comp[which(tech_comp$count_v1_1>=10),]
# rep_comp_cutoff<-rep_comp[which(rep_comp$count_v2_155>=10),]
# rnaseq_comp_cutoff<-rnaseq_comp[which(rnaseq_comp$count_rnaseq_155>=10),]

print(cor(tech_comp$count_v1_1,tech_comp$count_v2_1),method="spearman")
print(nrow(tech_comp))
print(cor(rep_comp$count_v2_155,rep_comp$count_v2_1),method="spearman")
print(nrow(rep_comp))

tech_comp<-tech_comp[sample(1:nrow(tech_comp),10000),]
rep_comp<-rep_comp[sample(1:nrow(rep_comp),10000),]

# p1<-ggplot(tech_comp,aes(x=count_v1_1,y=count_v2_1))+geom_point(size=0.5,alpha=0.1)+
# scale_x_continuous(name=NULL,trans="log10",limits=c(0.01,100000),breaks=log_ticks(0.01,100000)[[2]],labels=log_ticks_exp(0.01,100000)[[1]])+
# scale_y_continuous(name=NULL,trans="log10",limits=c(0.01,100000),breaks=log_ticks(0.01,100000)[[2]],labels=log_ticks_exp(0.01,100000)[[1]])+theme_tim()+geom_abline(linetype="dashed",color="grey")#+expand_limits(x = 0, y = 0)

# p2<-ggplot(rep_comp,aes(x=count_v2_155,y=count_v2_1))+geom_point(size=0.5,alpha=0.1)+
# scale_x_continuous(name=NULL,trans="log10",limits=c(.01,100000),breaks=log_ticks(.01,100000)[[2]],labels=log_ticks_exp(.01,100000)[[1]])+
# scale_y_continuous(name=NULL,trans="log10",limits=c(0.01,100000),breaks=log_ticks(0.01,100000)[[2]],labels=log_ticks_exp(0.01,100000)[[1]])+theme_tim()+
# geom_abline(linetype="dashed",color="grey")#+expand_limits(x = 0, y = 0)


p1<-ggplot(tech_comp,aes(x=log10(count_v1_1),y=log10(count_v2_1)))+geom_point(size=0.5,alpha=0.1)+
scale_x_continuous(name=NULL,limits=c(-2,5))+
scale_y_continuous(name=NULL,limits=c(-2,5))+theme_tim()+geom_abline(linetype="dashed",color="grey")#+expand_limits(x = 0, y = 0)

p2<-ggplot(rep_comp,aes(x=log10(count_v2_155),y=log10(count_v2_1)))+geom_point(size=0.5,alpha=0.1)+
scale_x_continuous(name=NULL,limits=c(-2,5))+
scale_y_continuous(name=NULL,limits=c(-2,5))+theme_tim()+
geom_abline(linetype="dashed",color="grey")#+expand_limits(x = 0, y = 0)

# p3<-ggplot(rnaseq_comp_cutoff,aes(x=count_rnaseq_155,y=count_rnaseq_1))+geom_point(size=0.5,alpha=0.5)+
# scale_x_continuous(name=NULL,trans="log10",limits=c(10,100000),breaks=breaks,labels=labels)+
# scale_y_continuous(name=NULL,trans="log10",limits=c(0.01,100000),breaks=breaks,labels=labels)+theme_tim()+
# geom_abline(linetype="dashed",color="grey")#+expand_limits(x = 0, y = 0)

ggsave(plot=p1,file="figures/version_6/FigS2Epanel1.png",width=2,height=2)
ggsave(plot=p2,file="figures/version_6/FigS2Epanel2.png",width=2,height=2)
# ggsave(plot=p3,file="SXX3.pdf",width=3,height=3)