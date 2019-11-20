library(plyr)
library(ggplot2)
library(reshape2)
library(cowplot)
source("/lab/solexa_bartel/teisen/RNAseq/Scripts/general/ggplot_theme.R")

log_ticks<-function(start,end){
  tic_num<-log10(end)-log10(start)+1
  labels<-rep(" ",10*tic_num)
  labels[seq(10,(tic_num*10),10)]<-10^(log10(start):log10(end))*10
  breaks<-c(matrix(1:10)%*%t(matrix(10^(log10(start):log10(end)))))
  return(list(labels,breaks))
}

halflife_binnings = TRUE

expr<-read.table("/lab/solexa_bartel/teisen/RNAseq/kinetics_of_translation/actD_RNAseq/analysis/minus_counts_10rpm_cutoff_std3_normalized_0hr_cutoff_only.txt",head=TRUE)
colnames(expr)<-c("accession","actD_0hr_expr","actD_1hr_expr","actD_3hr_expr","actD_7hr_expr","actD_15h_expr")
tails<-read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-1_actD_samples/analysis_20171030/model_input_files/means/background_subtracted_mean_tails_minus_50tags.txt",head=TRUE)
tails<-merge(expr,tails,by="accession")

# actD_0hr<-tails[,c(1,2,7)]
# actD_1hr<-tails[,c(1,3,8)]
# actD_3hr<-tails[,c(1,4,9)]
# actD_7hr<-tails[,c(1,5,10)]
# actD_15h<-tails[,c(1,6,11)]

print("Reading data files ...")
actD_0hr <- read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-1_actD_samples/analysis_20171030/model_input_files/global_normalized/actD_0hr_minus_background_subtracted_single_tag_v6.txt",sep="\t")
actD_1hr <- read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-1_actD_samples/analysis_20171030/model_input_files/global_normalized/actD_1hr_minus_background_subtracted_single_tag_v6.txt",sep="\t")
actD_3hr <- read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-1_actD_samples/analysis_20171030/model_input_files/global_normalized/actD_3hr_minus_background_subtracted_single_tag_v6.txt",sep="\t")
actD_7hr <- read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-1_actD_samples/analysis_20171030/model_input_files/global_normalized/actD_7hr_minus_background_subtracted_single_tag_v6.txt",sep="\t")
actD_15h <- read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-1_actD_samples/analysis_20171030/model_input_files/global_normalized/actD_15h_minus_background_subtracted_single_tag_v6.txt",sep="\t")
print("Calculating mean tails ...")
# actD_0hr_means <- apply(actD_0hr,1,function(x){weighted.mean(x=0:250,w=x[-c(252)])})
# actD_1hr_means <- apply(actD_1hr,1,function(x){weighted.mean(x=0:250,w=x[-c(252)])})
# actD_3hr_means <- apply(actD_3hr,1,function(x){weighted.mean(x=0:250,w=x[-c(252)])})
# actD_7hr_means <- apply(actD_7hr,1,function(x){weighted.mean(x=0:250,w=x[-c(252)])})
# actD_15h_means <- apply(actD_15h,1,function(x){weighted.mean(x=0:250,w=x[-c(252)])})

actD_0hr$means <- apply(actD_0hr,1,function(x){weighted.mean(x=8:250,w=x[-c(1:8,252)])})
actD_1hr$means <- apply(actD_1hr,1,function(x){weighted.mean(x=8:250,w=x[-c(1:8,252)])})
actD_3hr$means <- apply(actD_3hr,1,function(x){weighted.mean(x=8:250,w=x[-c(1:8,252)])})
actD_7hr$means <- apply(actD_7hr,1,function(x){weighted.mean(x=8:250,w=x[-c(1:8,252)])})
actD_15h$means <- apply(actD_15h,1,function(x){weighted.mean(x=8:250,w=x[-c(1:8,252)])})


actD_0hr<-actD_0hr[,c(252,253)]
actD_1hr<-actD_1hr[,c(252,253)]
actD_3hr<-actD_3hr[,c(252,253)]
actD_7hr<-actD_7hr[,c(252,253)]
actD_15h<-actD_15h[,c(252,253)]

print(nrow(actD_0hr))
print(nrow(actD_1hr))
print(nrow(actD_3hr))
print(nrow(actD_7hr))
print(nrow(actD_15h))

colnames(actD_0hr)<-c("accession","tail_length")
colnames(actD_1hr)<-c("accession","tail_length")
colnames(actD_3hr)<-c("accession","tail_length")
colnames(actD_7hr)<-c("accession","tail_length")
colnames(actD_15h)<-c("accession","tail_length")


actD_0hr<-merge(actD_0hr,expr[,c(1,2)])
actD_1hr<-merge(actD_1hr,expr[,c(1,3)])
actD_3hr<-merge(actD_3hr,expr[,c(1,4)])
actD_7hr<-merge(actD_7hr,expr[,c(1,5)])
actD_15h<-merge(actD_15h,expr[,c(1,6)])

colnames(actD_0hr)[3]<-"expr"
colnames(actD_1hr)[3]<-"expr"
colnames(actD_3hr)[3]<-"expr"
colnames(actD_7hr)[3]<-"expr"
colnames(actD_15h)[3]<-"expr"

actD_0hr$sample<-"h0"
actD_1hr$sample<-"h1"
actD_3hr$sample<-"h3"
actD_7hr$sample<-"h7"
actD_15h$sample<-"h15"

m<-rbind(
actD_0hr,
actD_1hr,
actD_3hr,
actD_7hr,
actD_15h)

#write.table(m,file="processed_files/minus_eluate_samples_50tag_cutoff_means.txt",quote=FALSE,sep="\t",row.names=FALSE)


#m<-m[which(m$sample %in% c("0","1","3")),]
print("Calculating gene bins ...")

ClosestNumber <- function(x,your.number){which(abs(x-your.number)==min(abs(x-your.number)))}

	
hl <- read.table("processed_files/halflifeComparisons/miR-155_minu_PAL_halflives_logspace_global_offset_SS.txt",head=TRUE,sep="\t")
hl <- hl[which(hl$accession %in% m$accession),]
hl$subset<-"residual"
hl<-hl[order(hl$halflife_t),]
hl<-hl[which(hl$accession %in% hl$accession),]

print(nrow(hl))

"hr1" <- ClosestNumber(hl$halflife,0.2427001)
"hr2" <- ClosestNumber(hl$halflife,0.9705460)
"hr3" <- ClosestNumber(hl$halflife,2.1707216)
"hr4" <- ClosestNumber(hl$halflife,4.3490184)
"hr5" <- ClosestNumber(hl$halflife,11.4832559)


"hr1" <-(hr1-49):(hr1+50)
"hr2" <-(hr2-49):(hr2+50)
"hr3" <-(hr3-49):(hr3+50)
"hr4" <-(hr4-49):(hr4+50)
"hr5" <-(hr5-49):(hr5+50)



hl[hr1,]$subset<-"hr1"
hl[hr2,]$subset<-"hr2"
hl[hr3,]$subset<-"hr3"
hl[hr4,]$subset<-"hr4"
hl[hr5,]$subset<-"hr5"


#long_hl_genes<-"NM_026147"
m<-merge(m,hl,by="accession")



print("Plotting ...")

# genes_to_remove<-m[m$tail_length>250,]$accession
#print(genes_to_remove)

# m<-m[which(!m$accession %in% genes_to_remove),]
print(length(unique(m[m$subset!="residual",]$accession)))
m$sample<-factor(m$sample,levels=c("h0","h1","h3","h7","h15"))
m$subset<-factor(m$subset,levels=c("hr5","hr4","hr3","hr2","hr1","residual"))


print(length(unique(m[which(m$subset=="hr1"),]$accession)))
print(length(unique(m[which(m$subset=="hr2"),]$accession)))
print(length(unique(m[which(m$subset=="hr3"),]$accession)))
print(length(unique(m[which(m$subset=="hr4"),]$accession)))
print(length(unique(m[which(m$subset=="hr5"),]$accession)))

g_legend<-function(a.gplot){ 
  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)} 

m<-m[order(m$sample),]
# cols<-c("#ccece6","#99d8c9","#66c2a4","#2ca25f","#006d2c")
# write.table(m[which(!m$subset == "residual"),],file="/lab/solexa_bartel/teisen/Tail-seq/miR-1_actD_samples/analysis_20171030/processed_files/fig5C_genes_v3.txt",row.names=FALSE,sep="\t",quote=FALSE)
cols<-c("#D53F50","#F47E57","#FBE28B","#B4D895","#3288BD")

levels(m$subset)=c(""," ","  ","   ","    ","     ")
p6<-ggplot(m[m$subset!="     ",],aes(x=tail_length,y=expr,color=sample))+
geom_path(alpha=.5,size=0.2,aes(group=accession))+
facet_wrap(~subset,ncol=6)+
scale_y_continuous(trans="log10",name=NULL,expand=c(0,0),limits=c(0.001,100),labels=log_ticks(0.001,1)[[1]],breaks=log_ticks(0.001,1)[[2]])+
scale_x_continuous(name="Tail Length (nt)",expand=c(0,0),limits=c(0,250))+
guides(colour = guide_legend(override.aes = list(alpha=1)))+
geom_point(size=0.1)+
scale_color_manual(values=cols)+
theme_tim()+
theme(axis.line.y=element_blank())+
geom_segment(aes(y=0.001,yend=10,x=0,xend=0),color="black",size=0.5)


data.segm<-data.frame(x=0,y=0.001,xend=0,yend=10,
                      subset="0.95")

levels(m$subset)=c("0.95","0.75","0.5","0.25","0.05","na")
p6a<-ggplot(m[m$subset!="na",],aes(x=tail_length,y=expr,color=sample))+
geom_path(alpha=.5,size=0.2,aes(group=accession))+
facet_wrap(~subset,ncol=6)+
scale_y_continuous(trans="log10",name=NULL,expand=c(0,0),limits=c(0.001,100),labels=log_ticks(0.001,1)[[1]],breaks=log_ticks(0.001,1)[[2]])+
scale_x_continuous(name="Tail Length (nt)",expand=c(0,0),limits=c(0,250))+
guides(colour = guide_legend(override.aes = list(alpha=1)))+
geom_point(size=0.1)+
scale_color_manual(values=cols)+
theme_tim()+
theme(axis.line.y=element_blank())+
theme(axis.line.y=element_blank(),plot.margin=unit(c(0,10,0,0), unit = "mm"))+
geom_segment(data=data.segm,aes(y=y,yend=yend,x=x,xend=xend),inherit.aes=FALSE,color="black",size=0.5)


msmooth <- m
smBin <- NULL
for(i in unique(msmooth$sample)){
  for(j in unique(msmooth$subset)){
    smBin <- rbind(smBin,data.frame(
      tail_length = median(msmooth[msmooth$sample==i & msmooth$subset==j,]$tail_length,na.rm=TRUE),
      expr = median(msmooth[msmooth$sample==i & msmooth$subset==j,]$expr,na.rm=TRUE),
      sample = i,
      subset = j))
}
} 

p7<-ggplot(m[m$subset!="     ",],aes(x=tail_length,y=expr,color=sample))+
geom_path(alpha=.5,size=0.2,aes(group=accession))+
facet_wrap(~subset,ncol=6)+
scale_y_continuous(trans="log10",name=NULL,expand=c(0,0),limits=c(0.001,100),labels=log_ticks(0.001,1)[[1]],breaks=log_ticks(0.001,1)[[2]])+
scale_x_continuous(name="Tail Length (nt)",expand=c(0,0),limits=c(0,250))+
guides(colour = guide_legend(override.aes = list(alpha=1)))+
geom_point(size=0.1)+
geom_path(data=smBin[smBin$subset!="     ",],aes(x=tail_length,y=expr,group=1),color="black",size=1,lineend = "round")+
scale_color_manual(values=cols)+
theme_tim()+
theme(axis.line.y=element_blank())+
geom_segment(aes(y=0.001,yend=10,x=0,xend=0),color="black",size=0.5)



#model plots
t2t_l <- read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-1_actD_samples/analysis_20171030/processed_files/fig5C_data_v8_Bounded.txt",head=TRUE,sep="\t")
colnames(t2t_l)[5]<-"accession"
t2t_l<-t2t_l[which(t2t_l$rs>5E-8),]
t2t_l$hl_bin<-factor(t2t_l$hl_bin,levels=c("hr5","hr4","hr3","hr2","hr1"))
levels(t2t_l$hl_bin)<-c("0.95","0.75","0.5","0.25","0.05")

colfunc <- colorRampPalette(c("black", "white"))
colfunc1 <- colorRampPalette(c(cols[1],cols[2]))
colfunc2 <- colorRampPalette(c(cols[2],cols[3]))
colfunc3 <- colorRampPalette(c(cols[3],cols[4]))
colfunc4 <- colorRampPalette(c(cols[4],cols[5]))

coldata <- data.frame(
  time = 1:max(t2t_l$time),
  color = c(colfunc(59),colfunc1(60),colfunc2(120),colfunc3(4*60),colfunc4(100)),stringsAsFactors=FALSE)

t2t_l <- merge(t2t_l,coldata,by="time")
colnames(t2t_l)[6]<-"cola"
# p4<-ggplot(t2t_l[which(t2t_l$time>=60),],aes(x=tl,y=rs,color=time,group=accession))+geom_line(size=0.2)+scale_x_continuous(limits=c(0,250),expand=c(0,0))+scale_y_continuous(trans="log10",limits=c(1e-8,1e-4),labels=log_ticks(1e-8,1e-4)[[1]],breaks=log_ticks(1e-8,1e-4)[[2]],name=NULL,expand=c(0,0))+facet_wrap(~hl_bin,nrow=1)+scale_colour_distiller(palette="Spectral",direction = 1)+theme_tim()+theme(strip.text.x=element_text(size=6))
nvec <- t2t_l[which(t2t_l$time>=60),]$cola
names(nvec) <- t2t_l[which(t2t_l$time>=60),]$time
t2t_l$rs = t2t_l$rs*1e5 #normalize y axis scales.
p4<-ggplot(t2t_l[which(t2t_l$time>=60),],aes(x=tl,y=rs,group=accession))+geom_line(size=0.2,aes(color=factor(time)))+scale_x_continuous(limits=c(0,250),expand=c(0,0))+scale_y_continuous(trans="log10",limits=c(0.001,100),labels=log_ticks(0.001,1)[[1]],breaks=log_ticks(0.001,1)[[2]],name=NULL,expand=c(0,0))+facet_wrap(~hl_bin,nrow=1)+scale_color_manual(values=nvec)+theme(strip.text.x=element_text(size=6))+theme_tim()+theme(axis.line.y=element_blank())+geom_segment(aes(y=0.001,yend=10,x=0,xend=0),color="black",size=0.5)

levels(t2t_l$hl_bin)=c(""," ","  ","   ","    ")
data.segm2<-data.frame(x=0,y=0.001,xend=0,yend=10,
                      hl_bin="")
p4a<-ggplot(t2t_l[which(t2t_l$time>=60),],aes(x=tl,y=rs,group=accession))+
geom_line(size=0.2,aes(color=factor(time)))+
scale_x_continuous(limits=c(0,250),expand=c(0,0))+
scale_y_continuous(trans="log10",limits=c(0.001,100),labels=log_ticks(0.001,1)[[1]],breaks=log_ticks(0.001,1)[[2]],name=NULL,expand=c(0,0))+
facet_wrap(~hl_bin,nrow=1)+
scale_color_manual(values=nvec)+
theme(strip.text.x=element_text(size=6))+
theme_tim()+
theme(axis.line.y=element_blank(),plot.margin=unit(c(0,10,0,0), unit = "mm"))+
geom_segment(data=data.segm2,aes(y=y,yend=yend,x=x,xend=xend),inherit.aes=FALSE,color="black",size=0.5)


legend <- get_legend(ggplot(t2t_l[which(t2t_l$time>=60),],aes(x=tl,y=rs,group=accession))+geom_line(size=0.2,aes(color=factor(time)))+scale_x_continuous(limits=c(0,250),expand=c(0,0))+scale_y_continuous(trans="log10",limits=c(0.001,10),labels=log_ticks(0.001,10)[[1]],breaks=log_ticks(0.001,10)[[2]],name=NULL,expand=c(0,0))+facet_wrap(~hl_bin,nrow=1)+theme(strip.text.x=element_text(size=6))+scale_color_manual(values=nvec))

msmooth <- t2t_l[t2t_l$time %in% (round(c(0,1,3,7,15)*60/1.73)+59),] #only at times indicated in legend
smBinList <- NULL
counter = 0
for(i in unique(msmooth$hl_bin)){
  for(j in unique(msmooth$time)){
    counter = counter+1
    smBinList[[counter]] <- data.frame(
      tl = median(msmooth[msmooth$hl_bin==i & msmooth$time==j,]$tl,na.rm=TRUE),
      rs = median(msmooth[msmooth$hl_bin==i & msmooth$time==j,]$rs,na.rm=TRUE),
      hl_bin = i,
      time = j)
}
} 
smBinData = do.call(rbind, smBinList)

# counter = 0
# smBinDataList2 = NULL
# for(i in unique(smBinData$hl_bin)){
#   bin = 50
#   counter = counter +1
#   tempData = smBinData[smBinData$hl == i,]
#   tempData$hl_bin = as.numeric(as.character(tempData$hl_bin))
#   tempData = tempData[order(tempData$time),]
#   smBinDataList2[[counter]]=t(sapply(1:(nrow(tempData)-bin),function(x){apply(tempData[x:(x+bin),],2,function(x){median(x,na.rm=TRUE)})}))
# }
# smBinData2 = data.frame(do.call(rbind, smBinDataList2))
# smBinData2$hl_bin = as.factor(smBinData2$hl_bin)

p5<-ggplot(t2t_l[which(t2t_l$time>=60),],aes(x=tl,y=rs,group=accession))+
geom_line(size=0.2,aes(color=factor(time)))+
geom_path(data=smBinData,aes(x=tl,y=rs,group = 1),color="black",size=1,lineend = "round")+
scale_x_continuous(limits=c(0,250),expand=c(0,0))+
scale_y_continuous(trans="log10",limits=c(0.001,100),labels=log_ticks(0.001,1)[[1]],breaks=log_ticks(0.001,1)[[2]],name=NULL,expand=c(0,0))+
facet_wrap(~hl_bin,nrow=1)+
scale_color_manual(values=nvec)+
theme(strip.text.x=element_text(size=6))+
theme_tim()+
theme(axis.line.y=element_blank())+
geom_segment(aes(y=0.001,yend=10,x=0,xend=0),color="black",size=0.5)


# legend <- get_legend(ggplot(t2t_l[which(t2t_l$time>=60),],aes(x=tl,y=rs,color=time,group=accession))+geom_line(size=0.2)+scale_x_continuous(limits=c(0,250),expand=c(0,0))+scale_y_continuous(trans="log10",limits=c(1e-8,1e-4),labels=log_ticks(.001,10)[[1]],breaks=log_ticks(1e-8,1e-4)[[2]],name=NULL,expand=c(0,0))+facet_wrap(~hl_bin,nrow=1)+scale_colour_distiller(palette="Spectral",direction = 1))+scale_fill_continuous(breaks=c(60,180))

 #BREAK STATEMENT HERE
t2t_lrsq <- t2t_l[t2t_l$time %in% (round(c(0,1,3,7,15)*60/1.73)+59),]
t2t_lrsq$cola <- NULL
t2t_lrsq$hl_bin <- NULL
t2t_lrsq$sample = "none"
t2t_lrsq[t2t_lrsq$time==59,]$sample = "h0"
t2t_lrsq[t2t_lrsq$time==94,]$sample = "h1"
t2t_lrsq[t2t_lrsq$time==163,]$sample = "h3"
t2t_lrsq[t2t_lrsq$time==302,]$sample = "h7"
t2t_lrsq[t2t_lrsq$time==579,]$sample = "h15"
t2t_lrsq$time <- NULL

mcomp <- m
mcomp$alpha_t <- NULL
mcomp$beta_t <- NULL
mcomp$halflife_t <- NULL
mcomp$residual_t <- NULL
mcomp$offset <- NULL
mcomp$subset <- NULL
mboth <- merge(t2t_lrsq,mcomp,by=c("accession","sample"))
print(cor(mboth$expr,mboth$rs,method="spearman"))
print(cor(mboth$tl,mboth$tail_length,method="spearman"))
print(length(unique(mboth$accession)))

ggsave(plot=p4a,file="/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/figures/version_7/Fig5D.png",width=6.33,height=2.2)
ggsave(plot=p6a,file="/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/figures/version_7/Fig5E.png",width=6.33,height=2.2)

break
ggsave(plot=p5,file="/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/figures/version_6/Fig5D.png",width=6,height=2.2)
ggsave(plot=p7,file="/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/figures/version_6/Fig5E.png",width=6,height=2.2)

ggsave(plot=p4,file="/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/figures/version_6/Fig5D.pdf",width=6,height=2.2,useDingbats=FALSE)
ggsave(plot=p6,file="/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/figures/version_6/Fig5E.pdf",width=6,height=2.2,useDingbats=FALSE)
# ggsave(plot=legend,file="/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/figures/version_6/Fig5Dlegend.pdf",width=19,height=19,useDingbats=FALSE)


