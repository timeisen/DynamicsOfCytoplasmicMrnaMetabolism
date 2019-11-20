library(tidyverse)
library(scales)
source('/lab/solexa_bartel/teisen/RNAseq/Scripts/general/ggplot_theme.R')
source('/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/KinetTlnAnalysis/NorthernAnalysis/NorthernAnalysisHelpers.R')


eef2 <- read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/processed_files/northerns/eef2_lane_profile.txt",head=TRUE,sep = "\t",skip=5)
colnames(eef2)<-c("pixel",paste0("lane_",1:14))

# runAll <- function(){
# 	source("KinetTlnAnalysis/NorthernAnalysis/NorthernAnalysisThbs1.R")
# 	source("KinetTlnAnalysis/NorthernAnalysis/NorthernAnalysisActb.R")
# 	source("KinetTlnAnalysis/NorthernAnalysis/NorthernAnalysisEef2.R")
# }

p1 <- ggplot(eef2,aes(x=pixel,y=lane_1,labels=pixel))+geom_line()+theme_bw()
p2 <- ggplot(eef2,aes(x=pixel,y=lane_2,labels=pixel))+geom_line()+theme_bw()


CplusMapping = data.frame(
	nt = c(500,400,300,200,100),
	pixel = c(23,40,71,133,308)
	)

DecMapping = data.frame(
	nt = c(150,100,90,80,70,60,50,40,30),
	pixel = c(183,306,336,372,411,458,510,572,642)
	)

BothMarkers = data.frame(
	nt = c(500,400,300,200,150,100,90,80,70,60,50,40,30),
	pixel = c(23,40,71,133,183,307,336,372,411,458,510,572,642)
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

# binning <- function(tailsDF){
# 	tailsIntensity <- data.frame(rep(0,606),rep(0,606),rep("a",606),stringsAsFactors=FALSE)
# 	jCount = 0
# 	for(j in unique(tailsDF$sample)){
# 		for(i in 0:100){
# 			intensity = mean(tailsDF[which(tailsDF$tail_length>i*5 & tailsDF$tail_length <= i*5+5 & tailsDF$sample == j),]$intensity)
# 			tailsIntensity[i+1+jCount*101,1] <- i*5
# 			tailsIntensity[i+1+jCount*101,2] <- intensity
# 			tailsIntensity[i+1+jCount*101,3] <- as.character(j)
# 		}
# 		jCount = jCount+1
# 	}
# 	print(head(tailsIntensity))
# 	colnames(tailsIntensity) = c("tail_length","intensity","sample")
# 	tailsIntensity$sample <- as.factor(tailsIntensity$sample)
# 	tailsIntensity$sample<-factor(tailsIntensity$sample,levels=c("SS","8 h","4 h","2 h","1 h","40 min"))
# 	return(tailsIntensity)
# }


p3 <- ggplot(BothMarkers,aes(x=nt,y=log(pixel)))+geom_point()+geom_abline(slope=bm$coefficients[[2]],intercept=bm$coefficients[[1]])

eef2$nt <- conversion(eef2$pixel,Obj_loess)
gather_eef2 <- gather(eef2,key="lane",value="value",-pixel,-nt)

gather_eef2$bg <- rep(apply(eef2[which(eef2$nt>400),],2,mean)[2:15],each=668)


cols<-rev(c(
	"black",
	"#8856a7",
	"#2b8cbe",
	"#31a354",
	"#fec44f",
	"#f03b20"))

laneSample <- data.frame(lane = c("lane_3","lane_4","lane_5","lane_6","lane_7","lane_8"),
		sample = rev(c("SS","8 h","4 h","2 h","1 h","40 min")))

tails <- gather_eef2[gather_eef2$lane %in% c("lane_3","lane_4","lane_5","lane_6","lane_7","lane_8"),]
tails <- merge(tails,laneSample)
tails$tail_length = tails$nt - 59
tails$intensity = tails$value-tails$bg
tails$sample<-factor(tails$sample,levels=rev(c("SS","8 h","4 h","2 h","1 h","40 min")))
tails$tail_length<-as.numeric(as.character(tails$tail_length))
tails$rescale<-0
tails$intensity = tails$intensity/10
tails$rescale <- 1/tails$intensity

print(median(tails$intensity))
print(head(tails))
print(min(tails$intensity))
tails$offset<-0

tails[which(tails$sample=="SS"),]$offset <- tails[which(tails$sample=="SS"),]$rescale
tails[which(tails$sample=="8 h"),]$offset <- tails[which(tails$sample=="8 h"),]$rescale+500
tails[which(tails$sample=="4 h"),]$offset <- tails[which(tails$sample=="4 h"),]$rescale+400
tails[which(tails$sample=="2 h"),]$offset <- tails[which(tails$sample=="2 h"),]$rescale+300
tails[which(tails$sample=="1 h"),]$offset <- tails[which(tails$sample=="1 h"),]$rescale+200
tails[which(tails$sample=="40 min"),]$offset <- tails[which(tails$sample=="40 min"),]$rescale + 100

#tails$offset <-tails$rescale+100*(as.numeric(tails$sample)-1)
scalerange <- as.numeric(quantile(tails$rescale,p=c(0.05,.95),na.rm=TRUE))

print(scalerange)
gradientends <- scalerange + rep(c(0,5,4,3,2,1)*100, each=2)

cols_hm<-cols
colorends <- c("white", cols_hm[1], "white", cols_hm[2], "white", cols_hm[3],"white",cols_hm[4],"white",cols_hm[5],"white",cols_hm[6])
print("Plotting ...")

p13 <- DataProcess(read_tsv("figures/other/otherV9/NM_007907.txt"),as_tibble(tails))

ggsave(plot=p13,width=2,height=4,file="figures/other/otherV9/Eef2_lane_profile_20191002_normalized.pdf",useDingbats=FALSE)