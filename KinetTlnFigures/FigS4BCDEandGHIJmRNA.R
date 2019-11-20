#FigS4BCDEandGHIJ mRNA

library(tidyverse)
# library(ggplot2)
# options("scipen"=-3)
source("ggplot_theme.R")
#plotting parameters here
alpha = 0.2
size  = 0.5
shape = 16
mone <- read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/rate_constant_measurements/miR-155_minus_samples_UNLINKV3_V81H5_run3_high_precision_reformat.txt")
mtwo150<-read.table("rate_constant_measurements/miR-155_minus_samples_UNLINKV3_V81H5_run3_high_precision_TwoDea150.txt")
mtwo110<-read.table("rate_constant_measurements/miR-155_minus_samples_UNLINKV3_V81H5_run3_high_precision_TwoDea110.txt")


colnames(mone)<-c("accession","st_1","a_1","k_1","b_1","r_1")
colnames(mtwo110)<-c("accession","st_110","a_110","k1_110","k2_110","b_110","r_110")
colnames(mtwo150)<-c("accession","st_150","a_150","k1_150","k2_150","b_150","r_150")


mall <- merge(mone,mtwo110,by="accession")
mall <- merge(mall,mtwo150,by="accession")

annot<-read.table("/lab/solexa_bartel/teisen/RNAseq/Annotation_files/TAIL_seq/symbol_to_accession.txt",head=TRUE)
mall<-merge(mall,annot,by="accession")

print(max(mall$r_1/1E6))
print(max(mall$r_110/1E6))
print(max(mall$r_150/1E6))
print(min(mall$r_1/1E6))
print(min(mall$r_110/1E6))
print(min(mall$r_150/1E6))

# mall[which(mall$r_1/1E6   < 10),]$r_1 = 10*1E6
# mall[which(mall$r_110/1E6 < 10),]$r_110 = 10*1E6
# mall[which(mall$r_150/1E6 < 10),]$r_150 = 10*1E6

mall[which(mall$r_1/1E6   > 1E4),]$r_1 =   1E4*1E6
mall[which(mall$r_110/1E6 > 1E4),]$r_110 = 1E4*1E6
mall[which(mall$r_150/1E6 > 1E4),]$r_150 = 1E4*1E6


p9_110 <- ggplot(mall,aes(x=r_1/1E6,y=r_110/1E6))+geom_point(alpha = alpha,size = size, shape = shape)+scale_x_continuous(trans="log10",
	expand=c(0,0),
	limits=c(1E1,1E4),
	labels=fancy_scientific,
	breaks=log_ticks(1E1,1E4)[[2]])+
	scale_y_continuous(trans="log10",name=NULL,
	expand=c(0,0),
	limits=c(1E1,1E4),
	labels=fancy_scientific,
	breaks=log_ticks(1E1,1E4)[[2]])+
	geom_abline(slope=1,intercept=0,linetype="dashed",color="grey")+theme_tim()

p9_150 <- ggplot(mall,aes(x=r_1/1E6,y=r_150/1E6))+geom_point(alpha = alpha,size = size, shape = shape)+scale_x_continuous(trans="log10",
	expand=c(0,0),
	limits=c(1E1,1E4),
	labels=fancy_scientific,
	breaks=log_ticks(1E1,1E4)[[2]])+
	scale_y_continuous(trans="log10",name=NULL,
	expand=c(0,0),
	limits=c(1E1,1E4),
	labels=fancy_scientific,
	breaks=log_ticks(1E1,1E4)[[2]])+
	geom_abline(slope=1,intercept=0,linetype="dashed",color="grey")+theme_tim()
print("max")
print(max(mall$k1_110))
print(max(mall$k1_150))
print(max(mall$k2_110)) #
print(max(mall$k2_150)) #
print("min")
print(min(mall$k1_110)) #
print(min(mall$k1_150)) #
print(min(mall$k2_110))
print(min(mall$k2_150))

# mall[which(mall$k1_110>1000),]$k1_110 = 1000
mall[which(mall$k2_110>1000),]$k2_110 = 1000
# mall[which(mall$k1_110<0.01),]$k1_110 = 0.01
# mall[which(mall$k1_150<0.01),]$k1_150 = 0.01

# mall[which(mall$k1_150>1000),]$k1_150 = 1000
mall[which(mall$k2_150>1000),]$k2_150 = 1000
# mall[which(mall$k2_150>1000),]$k2_150 = 1000
# mall[which(mall$k2_150<0.01),]$k2_150 = 0.01
# mall[which(mall$k1_150<0.01),]$k1_150 = 0.01

p10_110 <- ggplot(mall,aes(x=k1_110,y=k2_110,lables=symbol))+
	geom_point(alpha = alpha,size = size, shape = shape) +
	scale_x_continuous(
		expand=c(0,0),
		trans="log10",
		limits=c(0.01,1000),
		labels=fancy_scientific,
		breaks=log_ticks(.01,1000)[[2]])+scale_y_continuous(
		expand=c(0,0),
		trans="log10",name=NULL,
		limits=c(0.01,1000),
		labels=fancy_scientific,
		breaks=log_ticks(.01,1000)[[2]])+
	geom_abline(slope=1,intercept=0,linetype="dashed",color="grey")+theme_tim()

p10_150 <- ggplot(mall,aes(x=k1_150,y=k2_150,lables=symbol)) + 
	geom_point(alpha = alpha,size = size, shape = shape)+
	scale_x_continuous(
		expand=c(0,0),
		trans="log10",
		limits=c(0.01,1000),
		labels=fancy_scientific,
		breaks=log_ticks(.01,1000)[[2]])+scale_y_continuous(
		expand=c(0,0),
		trans="log10",name=NULL,
		limits=c(0.01,1000),
		labels=fancy_scientific,
		breaks=log_ticks(.01,1000)[[2]])+
	geom_abline(slope=1,intercept=0,linetype="dashed",color="grey")+theme_tim()
print(nrow(mall))
print(cor(mall$k1_110,mall$k2_110,method="spearman"))
print(cor(mall$r_1,mall$r_110,method="spearman"))
print(cor(mall$k1_150,mall$k2_150,method="spearman"))
print(cor(mall$r_1,mall$r_150,method="spearman"))
ggsave(plot=p9_110,file="/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/figures/version_9/FigS7B_110.pdf",width=2,height=2)
ggsave(plot=p10_110,file="/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/figures/version_9/FigS7C_110.pdf",width=2,height=2)
ggsave(plot=p9_150,file="/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/figures/version_9/FigS7B_150.pdf",width=2,height=2)
ggsave(plot=p10_150,file="/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/figures/version_9/FigS7C_150.pdf",width=2,height=2)


# print(nrow(mall))
# print(cor(mall$r_1,mall$r_2,method="spearman"))
# print(cor(mall$k1_2,mall$k2_2,method="spearman"))

#######CODE FOR 7BC ends here

m1 <- read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/rate_constant_measurements/miR-1_minus_samples_UNLINKV3_V81H5_run3_high_precision_reformat.txt") #updated 2019 08 13.

# m155 <- read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/rate_constant_measurements/miR-155_minus_samples_UNLINKV75_run6_final_global_param.txt")
m155 <- read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/rate_constant_measurements/miR-155_minus_samples_UNLINKV3_V81H5_run3_high_precision_reformat.txt")
# m1552 <- read.table("rate_constant_measurements/miR-155_minus_samples_UNLINKV76_run3_last7ntRm_HYBRID.txt")
hl <- read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/processed_files/halflifeComparisons/miR-155_minu_PAL_halflives_logspace_global_offset_SS_RNAseqScaling.txt",head=TRUE,sep="\t")

colnames(m1)<-c("accession","st_1","a_1","k_1","b_1","residual_1")
colnames(m155)<-c("accession","st_155","a_155","k_155","b_155","residual_155")
# colnames(m1552)<-c("accession","st_155h","a_155h","k_155h","b_155h","residual_155h")

mhl <- merge(m155,hl,by="accession")
# mhl2 <- merge(m155,m1552,by="accession")
m<-merge(m1,m155)
m<-m[complete.cases(m),]
print(cor(m[,-1],method="spearman"))
p1<-ggplot(m,aes(x=st_1,y=st_155))+geom_point(size = size, shape = shape,alpha = alpha)+
scale_x_continuous(name=NULL,limits=c(0,250),expand=c(0,1), breaks = c(0,100, 150, 250))+
scale_y_continuous(name=NULL,limits=c(0,250),expand=c(0,1), breaks = c(0,100, 150, 250))+
theme_tim()+
geom_abline(linetype="dashed",color="grey")+expand_limits(x = 0, y = 0)


# mA <- m
if (length(m[which(m$a_1 > 1E-5),]$a_1)>0){m[which(m$a_1 > 1E-5),]$a_1 = 1E-5}
if (length(m[which(m$a_155 > 1E-5),]$a_155)>0){m[which(m$a_155 > 1E-5),]$a_155 = 1E-5}
if (length(m[which(m$a_1 < 1E-8),]$a_1)>0){m[which(m$a_1 < 1E-8),]$a_1 = 1E-8}
if (length(m[which(m$a_155 < 1E-8),]$a_155)>0){m[which(m$a_155 < 1E-8),]$a_155 = 1E-8}

p2<-ggplot(m,aes(x=a_1,y=a_155))+geom_point(size = size, shape = shape,alpha = alpha)+
scale_x_continuous(expand=c(0,0),name=NULL,trans="log10",limits=c(1E-8,1E-5),labels=fancy_scientific,breaks=log_ticks(1E-9,1E-4)[[2]])+
scale_y_continuous(expand=c(0,0),name=NULL,trans="log10",limits=c(1E-8,1E-5),labels=fancy_scientific,breaks=log_ticks(1E-9,1E-4)[[2]])+
theme_tim()+
geom_abline(linetype="dashed",color="grey")+expand_limits(x = 0, y = 0)

m$kPRIME_155 = m$k_155
if (length(m[which(m$k_1 > 1E2),]$k_1)>0){m[which(m$k_1 > 1E2),]$k_1 = 1E2}
if (length(m[which(m$kPRIME_155 > 1E2),]$kPRIME_155)>0){m[which(m$kPRIME_155 > 1E2),]$kPRIME_155 = 1E2}
if (length(m[which(m$k_1 < 1E-2),]$k_1)>0){m[which(m$k_1 < 1E-2),]$k_1 = 1E-2}


p3<-ggplot(m,aes(x=k_1,y=kPRIME_155))+geom_point(size = size, shape = shape,alpha = alpha)+
scale_x_continuous(name=NULL,trans="log10",limits=c(0.03,30),labels=fancy_scientific,breaks=log_ticks(1E-2,1E+2)[[2]],expand=c(0,0))+
scale_y_continuous(name=NULL,trans="log10",limits=c(0.03,30),labels=fancy_scientific,breaks=log_ticks(1E-2,1E+2)[[2]],expand=c(0,0))+
theme_tim()+
geom_abline(linetype="dashed",color="grey")

pLogEl155 = plogis(230, loc = 263.95156,scale = 11.05133)
pLogEl1 = plogis(230, loc = 265.44891,scale = 13.97311)


if (length(m[which(m$b_1*pLogEl1 > 3),]$b_1)>0){m[which(m$b_1*pLogEl1 > 3),]$b_1 = 3/pLogEl1}
if (length(m[which(m$b_155*pLogEl155 > 3),]$b_155)>0){m[which(m$b_155*pLogEl155 > 3),]$b_155 = 3/pLogEl155}
if (length(m[which(m$b_1*pLogEl1 < 0.003),]$b_1)>0){m[which(m$b_1*pLogEl1 < 0.003),]$b_1 = 0.003/pLogEl1}
if (length(m[which(m$b_155*pLogEl155 < 0.003),]$b_155)>0){m[which(m$b_155*pLogEl155 < 0.003),]$b_155 = 0.003/pLogEl155}

p3b<-ggplot(m,aes(x=b_1*pLogEl1,y=b_155*pLogEl155))+geom_point(size = size, shape = shape,alpha = alpha)+
scale_x_continuous(name=NULL,trans="log10",limits=c(0.003,3),labels=fancy_scientific,breaks=log_ticks(1E-5,1E+2)[[2]],expand=c(0,0))+
scale_y_continuous(name=NULL,trans="log10",limits=c(0.003,3),labels=fancy_scientific,breaks=log_ticks(1E-5,1E+2)[[2]],expand=c(0,0))+
theme_tim()+
geom_abline(linetype="dashed",color="grey")

# mhl[which(mhl$alpha_t < 1E-8),]$alpha_t = 1E-8
mhl[which(mhl$a_155 < 1E-8),]$a_155 = 1E-8
mhl[which(mhl$alpha_t > 1E-5),]$alpha_t = 1E-5
mhl[which(mhl$a_155 > 1E-5),]$a_155 = 1E-5

p4<-ggplot(mhl,aes(x=alpha_t,y=a_155))+geom_point(size = size, shape = shape,alpha = alpha)+
scale_x_continuous(name=NULL,trans="log10",limits=c(1E-8,1E-5),labels=fancy_scientific,breaks=log_ticks(1E-8,1E-5)[[2]],expand=c(0,0))+
scale_y_continuous(name=NULL,trans="log10",limits=c(1E-8,1E-5),labels=fancy_scientific,breaks=log_ticks(1E-8,1E-5)[[2]],expand=c(0,0))+
theme_tim()+geom_abline(linetype="dashed",color="grey")

print("replicate correlation")
print(nrow(m))
print(cor(m[,-1],method="spearman"))

print("alpha")
print(nrow(mhl))
print(cor(mhl$alpha,mhl$a_155,method="spearman"))
print("hl deadenylation cor")
print(cor(mhl$halflife_t,mhl$k_155,method = "spearman"))
print(nrow(mhl))

bootstrap <- read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/processed_files/bootstrap/all_rates_means_sd_V81_H5_run3.txt",head=TRUE)
bootstrap$rate <- factor(bootstrap$rate,levels=c("st","a","k","b"))
# p5.1 <- ggplot(bootstrap[which(bootstrap$rate == "st"),],aes(x=sd/mean))+stat_ecdf()+facet_wrap(~rate)+scale_x_continuous(trans="log10",limits=c(1E-6,1E-1),labels=log_ticks(1E-6,1E-1)[[1]],breaks=log_ticks(1E-6,1E-1)[[2]])+theme_tim()+scale_y_continuous(name=NULL,expand=c(0,0),limits=c(0,1))
# p5.2 <- ggplot(bootstrap[which(bootstrap$rate == "a"),],aes(x=sd/mean))+stat_ecdf()+facet_wrap(~rate)+scale_x_continuous(trans="log10",limits=c(1E-6,1E-1),labels=log_ticks(1E-6,1E-1)[[1]],breaks=log_ticks(1E-6,1E-1)[[2]])+theme_tim()+scale_y_continuous(name=NULL,expand=c(0,0),limits=c(0,1))
# p5.3 <- ggplot(bootstrap[which(bootstrap$rate == "k"),],aes(x=sd/mean))+stat_ecdf()+facet_wrap(~rate)+scale_x_continuous(trans="log10",limits=c(1E-6,1E-1),labels=log_ticks(1E-6,1E-1)[[1]],breaks=log_ticks(1E-6,1E-1)[[2]])+theme_tim()+scale_y_continuous(name=NULL,expand=c(0,0),limits=c(0,1))

p5all <- ggplot(bootstrap,aes(x=sdlg,color=rate))+stat_ecdf()+scale_x_continuous(trans="log10",labels=fancy_scientific,breaks=log_ticks(1E-3,100)[[2]],expand=c(0,0))+theme_tim()+scale_y_continuous(name=NULL,expand=c(0,0),limits=c(0,1))+coord_cartesian(xlim=c(1E-3,1))+theme(axis.text.x=element_text(angle=45, hjust=1))
print("quantiles for bootstrap analysis")
print("st:")
print(quantile(bootstrap[which(bootstrap$rate == "st"),]$sdlg,probs = c(0.1,0.9)))
print("a:")
print(quantile(bootstrap[which(bootstrap$rate == "a"),]$sdlg,probs = c(0.1,0.9)))
print("k:")
print(quantile(bootstrap[which(bootstrap$rate == "k"),]$sdlg,probs = c(0.1,0.9)))
print("b:")
print(quantile(bootstrap[which(bootstrap$rate == "b"),]$sdlg,probs = c(0.1,0.9)))



randI <- read.table("rate_constant_measurements/miR-155_minus_samples_UNLINKV3_V81H5_run3_high_precision_RandI_reformat.txt",head=TRUE)
# randomized_param <- data.frame(randI$accession,randI$st_sd/randI$st_mean,randI$a_sd/randI$a_mean,randI$k_sd/randI$k_mean,randI$b_sd/randI$b_mean)
randomized_param <- data.frame(randI$accession,randI$st_sd,randI$a_sd,randI$k_sd,randI$b_sd)
colnames(randomized_param) <- c("accession","st","a","k","b")
rp <- gather(randomized_param,rate,cv,-accession)
rp$rate <- factor(rp$rate,levels=c("st","a","k","b"))

print("quantiles for randI analysis")
print("st:")
quantile(randomized_param$st,probs = c(0.1,0.9))
print("a:")
quantile(randomized_param$a,probs = c(0.1,0.9))
print("k:")
quantile(randomized_param$k,probs = c(0.1,0.9))
print("b:")
quantile(randomized_param$b,probs = c(0.1,0.9))

p6.1 <- ggplot(rp[which(rp$rate=="st"),],aes(x=cv))+stat_ecdf()+facet_wrap(~rate)+scale_x_continuous(trans="log10",limits=c(1E-6,1E-1),labels=log_ticks(1E-6,1E-1)[[1]],breaks=log_ticks(1E-6,1E-1)[[2]])+theme_tim()+scale_y_continuous(name=NULL,expand=c(0,0),limits=c(0,1))
p6.2 <- ggplot(rp[which(rp$rate=="a"),],aes(x=cv))+stat_ecdf()+facet_wrap(~rate)+scale_x_continuous(trans="log10",limits=c(1E-6,1E-1),labels=log_ticks(1E-6,1E-1)[[1]],breaks=log_ticks(1E-6,1E-1)[[2]])+theme_tim()+scale_y_continuous(name=NULL,expand=c(0,0),limits=c(0,1))
p6.3 <- ggplot(rp[which(rp$rate=="k"),],aes(x=cv))+stat_ecdf()+facet_wrap(~rate)+scale_x_continuous(trans="log10",limits=c(1E-6,1E-1),labels=log_ticks(1E-6,1E-1)[[1]],breaks=log_ticks(1E-6,1E-1)[[2]])+theme_tim()+scale_y_continuous(name=NULL,expand=c(0,0),limits=c(0,1))

p6all <- ggplot(rp,aes(x=cv,color=rate))+stat_ecdf()+scale_x_continuous(trans="log10",labels=fancy_scientific,breaks=log_ticks(1E-6,.1)[[2]],expand=c(0,0))+theme_tim()+scale_y_continuous(name=NULL,expand=c(0,0),limits=c(0,1))+coord_cartesian(xlim=c(1E-6,1E-3))+theme(axis.text.x=element_text(angle=45, hjust=1))

# if (length(m[which(m$k_155 < 1E-1),]$k_155)>0){m[which(m$k_155 < 1E-1),]$k_155 = 1E-1}
# m$shape = "\u25D7"

# m[which(m$b_155*pLogEl > 1E2),]$shape = "\u25D6"
# m[which(m$b_155*pLogEl < 1E-3),]$shape = "\u25D7"
# m[which(m$b_155*pLogEl > 1E2),]$b_155 = 1E2/pLogEl
m155[which(m155$b_155*pLogEl155 < 0.003),]$b_155 = 0.003/pLogEl155
m155[which(m155$b_155*pLogEl155 > 3),]$b_155 = 3/pLogEl155

m155[which(m155$k_155 > 30),]$k_155 = 30
m155[which(m155$k_155 < 0.03),]$k_155 = 0.03

p7<-ggplot(m,aes(x=b_155,y=k_155))+geom_point(alpha = alpha, shape = shape, size = size)+
scale_x_continuous(name=NULL,trans="log10",limits=c(1E-1,1E+3),labels=fancy_scientific,breaks=log_ticks(1E-1,1E+3)[[2]],expand=c(0,0))+
scale_y_continuous(name=NULL,trans="log10",limits=c(1E-1,1E+3),labels=fancy_scientific,breaks=log_ticks(1E-1,1E+3)[[2]],expand=c(0,0))+
theme_tim()#+theme(axis.text.x=element_text(angle=45, hjust=1))
print("plotting 4F")

p7x4<-ggplot(m155,aes(x=b_155*pLogEl155,y=k_155))+geom_point(size = size,alpha = alpha, shape  = shape)+
scale_x_continuous(name=NULL,trans="log10",limits=c(0.003,3),labels=log_ticks(1E-3,1E+2)[[1]],breaks=log_ticks(1E-3,1E+2)[[2]],expand=c(0,0))+
scale_y_continuous(name=NULL,trans="log10",limits=c(0.03,30),labels=log_ticks(1E-2,1E+3)[[1]],breaks=log_ticks(1E-2,1E+3)[[2]],expand=c(0,0))+
theme_tim()+
geom_abline(linetype="dashed",color="grey")
ggsave(plot=p7x4,file="figures/version_9/Fig4F.pdf",width=2,height=2,useDingbats=FALSE)
print("dea dcp cor")
print(cor(m155$b_155*pLogEl155,m155$k_155*pLogEl155,method = 'spearman'))

# mhl2 <- mhl2[,c(1,4,5,9,10)]
# mhl2df1 <- data.frame(mhl2$accession,mhl2$k_155,mhl2$b_155)
# mhl2df1$data <- "PALseq"
# mhl2df2 <- data.frame(mhl2$accession,mhl2$k_155h,mhl2$b_155h)
# mhl2df2$data <- "Hybrid"

# colnames(mhl2df1) <- c("accession","k","b","data")
# colnames(mhl2df2) <- c("accession","k","b","data")
# mhl2 <- rbind(mhl2df1,mhl2df2)

# p7x5<-ggplot(mhl2,aes(x=b*0.02016849,y=k))+geom_point(size=1,alpha=0.1,shape=16)+
# scale_x_continuous(name=NULL,trans="log10",limits=c(1E-3,1E+2),labels=fancy_scientific,breaks=log_ticks(1E-3,1E+2)[[2]],expand=c(0,0))+
# scale_y_continuous(name=NULL,trans="log10",limits=c(1E-1,1E+3),labels=fancy_scientific,breaks=log_ticks(1E-1,1E+3)[[2]],expand=c(0,0))+
# theme_tim()+
# facet_wrap(~data)+
# geom_abline(linetype="dashed",color="grey")
# # ggsave(plot=p7x5,file="figures/other/FigS7IrevisedHyrbid.pdf",width=4,height=2,useDingbats=FALSE)




if (length(mhl[which(mhl$halflife_t > 1E2),]$halflife_t)>0){mhl[which(mhl$halflife_t > 1E2),]$halflife_t = 1E2}
if (length(mhl[which(mhl$halflife_t < 1E-3),]$halflife_t)>0){mhl[which(mhl$halflife_t < 1E-3),]$halflife_t = 1E-3}
if (length(mhl[which(mhl$k_155 > 1E2),]$k_155)>0){mhl[which(mhl$k_155 > 1E2),]$k_155 = 1E2}
mhl <- merge(mhl,annot,by="accession")
options("scipen"=5)

mhl[mhl$halflife_t < 0.1,]$halflife_t = 0.1
p8 <- ggplot(mhl,aes(y=halflife_t,x=k_155,labels=symbol))+
geom_point(size = size, shape = shape,alpha = alpha)+
scale_y_continuous(name=NULL,trans="log10",limits=c(1E-1,1E+2),labels = log_ticks(1E-1,1E+3)[[1]],breaks=log_ticks(1E-1,1E+3)[[2]])+
scale_x_continuous(name=NULL,trans="log10",limits=c(0.01,1E+2),labels = log_ticks(1E-2,1E+2)[[1]],breaks=log_ticks(1E-2,1E+2)[[2]],expand=c(0,0)) +
geom_segment(aes(x = 0.01, xend = 0.01, y = 1E-1, yend = 1E2)) +
theme_tim() + 
theme(axis.line.y = element_blank()) 

p8a <- ggplot(mhl,aes(x=halflife_t,y=k_155,labels=symbol))+geom_point()+
scale_x_continuous(name=NULL,trans="log10",limits=c(1E-3,1E+2),labels=fancy_scientific,breaks=log_ticks(1E-3,1E+3)[[2]],expand=c(0,0))+
scale_y_continuous(name=NULL,trans="log10",limits=c(0.01,1E+2),labels=fancy_scientific,breaks=log_ticks(1E-2,1E+2)[[2]],expand=c(0,0))+
theme_tim()+theme(axis.text.x=element_text(angle=45, hjust=1))

ggsave(plot=p1,file="figures/version_9/FigS7Gpanel1.pdf",width=1.5,height=1.5,useDingbats=FALSE)
ggsave(plot=p2,file="figures/version_9/FigS7Gpanel2.pdf",width=1.5,height=1.5,useDingbats=FALSE)
ggsave(plot=p3,file="figures/version_9/FigS7Gpanel3.pdf",width=1.5,height=1.5,useDingbats=FALSE)
ggsave(plot=p3b,file="figures/version_9/FigS7Gpanel4.pdf",width=1.5,height=1.5,useDingbats=FALSE)
ggsave(plot=p4,file="figures/version_9/FigS7H.pdf",width=2,height=2,useDingbats=FALSE)
ggsave(plot=p5all,file="figures/version_9/FigS7F.pdf",width=2,height=2,useDingbats=FALSE)
ggsave(plot=p6all,file="figures/version_9/FigS7E.pdf",width=2,height=2,useDingbats=FALSE)
ggsave(plot=p7,file="figures/version_9/FigS7I.pdf",width=2,height=2,useDingbats=FALSE)
ggsave(plot=p8,file="figures/version_9/Fig4G.pdf",width=2,height=2,useDingbats=FALSE)


