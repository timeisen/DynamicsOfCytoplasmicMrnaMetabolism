library(plyr)
library(ggplot2)
library(data.table)
library(tidyverse)

m <- NULL
for(i in 1:10){
	mTemp <- read_tsv(paste0("rate_constant_measurements/bootstrap/miR-155_minus_miR-155_minus_samples_UNLINKV3_V81_H5_Run3_BOOTSTRAP_",i,".txt"),col_names = FALSE)
	colnames(mTemp)[c(1,6)] <- c("accession","residual")
	m[[i]] <- data.frame(mTemp[order(mTemp$residual),] %>% distinct(accession,.keep_all = TRUE)) #keep the lowest residual.
	colnames(m[[i]])<-c("accession",paste0("st_",i),paste0("a_",i),paste0("k_",i),paste0("b_",i),paste0("r_",i))
}

all <- join_all(m,by="accession",type="full")
# all <- all[complete.cases(all),]

st <- all[,c(1,2,1:9*5+2)]
a  <- all[,c(1,3,1:9*5+3)]
k  <- all[,c(1,4,1:9*5+4)]
b  <- all[,c(1,5,1:9*5+5)]
r  <- all[,c(1,6,1:9*5+6)]

#Can't contain more than 8 NAs
st <- st[which(rowSums(is.na(st)) < 8),]
a  <- a[which(rowSums(is.na(a)) < 8),]
k  <- k[which(rowSums(is.na(k)) < 8),]
b  <- b[which(rowSums(is.na(b)) < 8),]
r  <- r[which(rowSums(is.na(r)) < 8),]

m <- NULL
m <- list(st,a,k,b)
# break
for(i in 1:4){
	m[[i]]$mean <-  apply(m[[i]][,-1],1,function(x){mean(x,na.rm=TRUE)})
	m[[i]]$sd <- apply(m[[i]][,-1],1,function(x){sd(x,na.rm=TRUE)})
	m[[i]] <- m[[i]][,c(1,12,13)]
}

m[[1]]$rate <- "st"
m[[2]]$rate <- "a"
m[[3]]$rate <- "k"
m[[4]]$rate <- "b"
all_met <- rbindlist(m)

# p1 <- ggplot(all_met,aes(x=sd/mean))+stat_ecdf()+facet_wrap(~rate)+scale_x_continuous(trans="log10")
write.table(all_met,file="processed_files/bootstrap/all_rates_means_sd_V81_H5_run3.txt",quote=FALSE,sep="\t",row.names=FALSE)