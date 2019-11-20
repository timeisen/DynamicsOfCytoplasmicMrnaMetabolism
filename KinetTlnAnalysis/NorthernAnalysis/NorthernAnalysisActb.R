library(tidyverse)
source('/lab/solexa_bartel/teisen/RNAseq/Scripts/general/ggplot_theme.R')
source('/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/KinetTlnAnalysis/NorthernAnalysis/NorthernAnalysisHelpers.R')

library(scales)

actb <- read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/processed_files/northerns/actb_lane_profile.txt",head=TRUE,sep = "\t",skip=5)
colnames(actb)<-c("pixel",paste0("lane_",1:14))


p1 <- ggplot(actb,aes(x=pixel,y=lane_1,labels=pixel))+geom_line()+theme_bw()
p2 <- ggplot(actb,aes(x=pixel,y=lane_2,labels=pixel))+geom_line()+theme_bw()


CplusMapping = data.frame(
	nt = c(500,400,300,200,100),
	pixel = c(65,81,113,176,351)
	)

DecMapping = data.frame(
	nt = c(150,100,90,80,70,60,50,40,30),
	pixel = c(224,348,378,413,455,501,552,612,686)
	)

BothMarkers = data.frame(
	nt = c(500,400,300,200,150,100,90,80,70,60,50,40,30),
	pixel = c(65,81,113,176,224,349.5,378,413,455,501,552,612,686)
	)

bm_dec <- lm(CplusMapping$nt~log(CplusMapping$pixel))
bm_cp <- lm(DecMapping$nt~log(DecMapping$pixel))
bm <- lm(BothMarkers$nt~log(BothMarkers$pixel))

p4 <- ggplot(DecMapping,aes(y=nt,x=log(pixel)))+geom_point()+geom_abline(slope=bm_cp$coefficients[[2]],intercept=bm_cp$coefficients[[1]])
p5 <- ggplot(CplusMapping,aes(y=nt,x=log(pixel)))+geom_point()+geom_abline(slope=bm_cp$coefficients[[2]],intercept=bm_cp$coefficients[[1]])
p6 <- ggplot(BothMarkers,aes(y=nt,x=log(pixel)))+geom_point()+geom_smooth()
Obj_loess <- loess(BothMarkers$nt~log(BothMarkers$pixel))


conversion <- function(pixel,Obj_loess){
	return(predict(Obj_loess,log(pixel)))
}


p3 <- ggplot(BothMarkers,aes(x=nt,y=log(pixel)))+geom_point()+geom_abline(slope=bm$coefficients[[2]],intercept=bm$coefficients[[1]])

actb$nt <- conversion(actb$pixel,Obj_loess)
gather_actb <- gather(actb,key="lane",value="value",-pixel,-nt)

p7 <- ggplot(actb[-1,],aes(x=nt,y=lane_1,labels=pixel))+geom_line()+theme_bw()+geom_vline(xintercept=c(100,200,300,400),linetype="dashed",color="grey")
p8 <- ggplot(actb[-1,],aes(x=nt,y=lane_2,labels=pixel))+geom_line()+theme_bw()+geom_vline(xintercept=c(150,100,90,80,70,60,50,40,30),linetype="dashed",color="grey")+scale_x_continuous(limits=c(0,400))

#83 or 62 for the two beaks in the +dT lanes



gather_actb$expr <- rep(c(
	0,
	0,
	142.68/1356.56,
	145.14/1059.97,
	63.97/218.07,
	86.32/99.95,
	99.81/57.97,
	135.54/50.94,
	142.68/1356.56,
	145.14/1059.97,
	63.97/218.07,
	86.32/99.95,
	99.81/57.97,
	135.54/50.94),each=728)

#apply(actb[which(actb$nt>400),],2,mean)
gather_actb$bg <- rep(c(10296.0921,9637.1988,9954.6948,10023.5882,9820.6865,9963.3971,9696.6056,9156.4046,9189.9659,9008.4757,8370.2820,8010.5402,7573.7801,7130.6537),each=728)



p9 <- ggplot(gather_actb[gather_actb$lane %in% c("lane_9","lane_10","lane_11","lane_12","lane_13","lane_14"),],aes(x=nt,y=value,color=lane,labels=nt))+geom_line()+theme_bw()

p10 <- ggplot(gather_actb[gather_actb$lane %in% c("lane_3","lane_4","lane_5","lane_6","lane_7","lane_8"),],aes(x=nt-72.5,y=(value-bg)*expr,color=lane))+geom_line()+theme_bw()

cols<-rev(c(
	"black",
	"#8856a7",
	"#2b8cbe",
	"#31a354",
	"#fec44f",
	"#f03b20"))

laneSample <- data.frame(lane = c("lane_3","lane_4","lane_5","lane_6","lane_7","lane_8"),
		sample = rev(c("SS","8 h","4 h","2 h","1 h","40 min")))

tails <- gather_actb[gather_actb$lane %in% c("lane_3","lane_4","lane_5","lane_6","lane_7","lane_8"),]
tails <- merge(tails,laneSample)
tails$tail_length = tails$nt - 72.5
tails$intensity = tails$value-tails$bg

tails$sample<-factor(tails$sample,levels=rev(c("SS","8 h","4 h","2 h","1 h","40 min")))

print("Plotting ...")


p13 <- DataProcess(read_tsv("figures/other/otherV9/NM_007393.txt"),as_tibble(tails))


# p11c <- plot_grid(p11b,p11a,align='v',labels=NULL,nrow=2,rel_heights=c(3,1))
ggsave(plot=p13,width=2,height=4,file="figures/other/otherV9/actb_lane_profile_20191002_normalized.pdf",useDingbats=FALSE)

