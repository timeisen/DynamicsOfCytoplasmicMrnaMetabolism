library(tidyverse)
library(scales)
source('/lab/solexa_bartel/teisen/RNAseq/Scripts/general/ggplot_theme.R')
source('/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/KinetTlnAnalysis/NorthernAnalysis/NorthernAnalysisHelpers.R')

thbs1 <- read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/processed_files/northerns/thbs1_lane_profile_20180712.txt",head=TRUE,sep = "\t",skip=5)
colnames(thbs1)<-c("pixel",paste0("lane_",1:14))


p1 <- ggplot(thbs1,aes(x=pixel,y=lane_1,labels=pixel))+geom_line()+theme_bw()
p2 <- ggplot(thbs1,aes(x=pixel,y=lane_2,labels=pixel))+geom_line()+theme_bw()


CplusMapping = data.frame(
	nt = c(400,300,200,100),
	pixel = c(67,98,162,331)
	)

DecMapping = data.frame(
	nt = c(150,100,90,80,70,60,50,40,30),
	pixel = c(207,329,359,392,428,470,514,571,632)
	)

BothMarkers = data.frame(
	nt = c(400,300,200,150,100,90,80,70,60,50,40,30),
	pixel = c(67,98,162,207,330,359,392,428,470,514,571,632)
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

thbs1$nt <- conversion(thbs1$pixel,Obj_loess)
gather_thbs1 <- gather(thbs1,key="lane",value="value",-pixel,-nt)

p7 <- ggplot(thbs1[-1,],aes(x=nt,y=lane_1,labels=pixel))+geom_line()+theme_bw()+geom_vline(xintercept=c(100,200,300,400),linetype="dashed",color="grey")
p8 <- ggplot(thbs1[-1,],aes(x=nt,y=lane_2,labels=pixel))+geom_line()+theme_bw()+geom_vline(xintercept=c(150,100,90,80,70,60,50,40,30),linetype="dashed",color="grey")+scale_x_continuous(limits=c(0,400))
p9 <- ggplot(gather_thbs1[gather_thbs1$lane %in% c("lane_3","lane_4","lane_5","lane_6","lane_7","lane_8"),],aes(x=nt-40,y=value))+geom_line()+theme_bw()+scale_x_continuous(limits=c(0,400),name="Tail Length (nt)")+facet_wrap(~lane,ncol=1)

# ggsave(plot=p9,width=4,height=10,file="figures/other/thbs1_lane_profile_20180712.pdf",useDingbats=FALSE)

##This needs to be updated from thbs1
# gather_thbs1$expr <- rep(c(
# 	0,
# 	0,
# 	142.68/1356.56,
# 	145.14/1059.97,
# 	63.97/218.07,
# 	86.32/99.95,
# 	99.81/57.97,
# 	135.54/50.94,
# 	142.68/1356.56,
# 	145.14/1059.97,
# 	63.97/218.07,
# 	86.32/99.95,
# 	99.81/57.97,
# 	135.54/50.94),each=728)


# gather_thbs1$bg <- rep(apply(thbs1[which(thbs1$nt>50 & thbs1$nt<100),],2,mean)[2:15],each=682)

gather_thbs1$bg <- rep(apply(thbs1[which(thbs1$nt > 350),],2,function(x){mean(x,na.rm = TRUE)})[2:15],each=682)
# gather_thbs1$bg <- rep(c(10296.0921,9637.1988,9954.6948,10023.5882,9820.6865,9963.3971,9696.6056,9156.4046,9189.9659,9008.4757,8370.2820,8010.5402,7573.7801,7130.6537),each=682)


cols<-rev(c(
	"black",
	"#8856a7",
	"#2b8cbe",
	"#31a354",
	"#fec44f",
	"#f03b20"))

laneSample <- data.frame(lane = c("lane_3","lane_4","lane_5","lane_6","lane_7","lane_8"),
		sample = rev(c("SS","8 h","4 h","2 h","1 h","40 min")))

tails <- gather_thbs1[gather_thbs1$lane %in% c("lane_3","lane_4","lane_5","lane_6","lane_7","lane_8"),]
tails <- merge(tails,laneSample)
tails$tail_length = tails$nt - 52 #The zero. 
tails$intensity = tails$value-tails$bg

p10 <- ggplot(tails, aes(x = tail_length, y = intensity)) + 
	geom_line() + 
	geom_vline(xintercept = 0, linetype = 'dashed', color = 'grey') + 
	geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey') +
	facet_wrap(~lane,ncol = 1) + 
	theme_tim_label()

tails$sample<-factor(tails$sample,levels=rev(c("SS","8 h","4 h","2 h","1 h","40 min")))

print("Plotting ...")
	
p13 <- DataProcess(read_tsv("figures/other/otherV9/NM_011580.txt"),as_tibble(tails))

	
ggsave(plot=p13,width=2,height=4,file="figures/other/otherV9/Thbs1_lane_profile_20181002_normalized.pdf",useDingbats=FALSE)