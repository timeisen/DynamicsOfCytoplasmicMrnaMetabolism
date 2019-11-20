#FigS2E.R mRNA
#u_analysis_metagene_miR
library(plyr)
library(ggplot2)
library(reshape2)
source("/lab/solexa_bartel/teisen/RNAseq/Scripts/general/ggplot_theme.R")

plot_miR<-function(miR){
	if(miR=="miR1"){
	miR1_NoSite_tpUTR<-read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/other_analyses/u_analysis/miR_analysis_2/nosite_m1_tpUTR.txt")
	miR1_NoSite_topsite<-read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/other_analyses/u_analysis/miR_analysis_2/nosite_m1_topsite.txt")	
	miR1_TopSite<-read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/miR_analysis/miR-1_3T3_25pct_down.txt")
	miR1_TpUtr<-read.table("/archive/bartel/2015.03.04-11263/solexa_bartel/eichhorn/Ribosome_profiling_paper/Sites/Mouse/miR-1/mmu-miR-1_tpUTR_oneormore.txt")

	# miR1_NoSite_tpUTR<-read.table("/archive/bartel/2015.03.04-11263/solexa_bartel/eichhorn/Ribosome_profiling_paper/Sites/Mouse/miR-1/mmu-miR-1_NoSite.txt")
	# miR1_NoSite_topsite <- read.table("/archive/bartel/2015.03.04-11263/solexa_bartel/eichhorn/Ribosome_profiling_paper/Sites/Mouse/miR-1/mmu-miR-1_NoSite.txt")
	
	TAGTGC<-read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/other_analyses/u_analysis/miR_analysis_3/TAGTGC-1_miR.txt")
	GCTACA<-read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/other_analyses/u_analysis/miR_analysis_3/GCTACA-1_miR.txt")
	AATCCG<-read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/other_analyses/u_analysis/miR_analysis_3/AATCCG-1_miR.txt")
	ACTGGT<-read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/other_analyses/u_analysis/miR_analysis_3/ACTGGT-1_miR.txt")
	TGTCAC<-read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/other_analyses/u_analysis/miR_analysis_3/TGTCAC-1_miR.txt")
	CGGTTA<-read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/other_analyses/u_analysis/miR_analysis_3/CGGTTA-1_miR.txt")
	GTCTAG<-read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/other_analyses/u_analysis/miR_analysis_3/GTCTAG-1_miR.txt")
	TGAACT<-read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/other_analyses/u_analysis/miR_analysis_3/TGAACT-1_miR.txt")
	CTAGTC<-read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/other_analyses/u_analysis/miR_analysis_3/CTAGTC-1_miR.txt")
	ATCGAA<-read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/other_analyses/u_analysis/miR_analysis_3/ATCGAA-1_miR.txt")
	GCAATT<-read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/other_analyses/u_analysis/miR_analysis_3/GCAATT-1_miR.txt")
	CAGCGT<-read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/other_analyses/u_analysis/miR_analysis_3/CAGCGT-1_miR.txt")}

	if(miR=="miR155"){
	miR1_NoSite_tpUTR<-read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/other_analyses/u_analysis/miR_analysis_2/nosite_m155_tpUTR.txt")
	miR1_NoSite_topsite<-read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/other_analyses/u_analysis/miR_analysis_2/nosite_m155_topsite.txt")	
	miR1_NoSite<-read.table("/archive/bartel/2015.03.04-11263/solexa_bartel/eichhorn/Ribosome_profiling_paper/Sites/Mouse/miR-155/mmu-miR-155_NoSite.txt")
	miR1_TopSite<-read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/miR_analysis/miR-155_3T3_25pct_down.txt")
	miR1_TpUtr<-read.table("/archive/bartel/2015.03.04-11263/solexa_bartel/eichhorn/Ribosome_profiling_paper/Sites/Mouse/miR-155/mmu-miR-155_tpUTR_oneormore.txt")

	# miR1_NoSite_tpUTR<-read.table("/archive/bartel/2015.03.04-11263/solexa_bartel/eichhorn/Ribosome_profiling_paper/Sites/Mouse/miR-155/mmu-miR-155_NoSite.txt")
	# miR1_NoSite_topsite <- read.table("/archive/bartel/2015.03.04-11263/solexa_bartel/eichhorn/Ribosome_profiling_paper/Sites/Mouse/miR-155/mmu-miR-155_NoSite.txt")

	TAGTGC<-read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/other_analyses/u_analysis/miR_analysis_3/TAGTGC-2_miR.txt")
	GCTACA<-read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/other_analyses/u_analysis/miR_analysis_3/GCTACA-2_miR.txt")
	AATCCG<-read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/other_analyses/u_analysis/miR_analysis_3/AATCCG-2_miR.txt")
	ACTGGT<-read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/other_analyses/u_analysis/miR_analysis_3/ACTGGT-2_miR.txt")
	TGTCAC<-read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/other_analyses/u_analysis/miR_analysis_3/TGTCAC-2_miR.txt")
	CGGTTA<-read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/other_analyses/u_analysis/miR_analysis_3/CGGTTA-2_miR.txt")
	GTCTAG<-read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/other_analyses/u_analysis/miR_analysis_3/GTCTAG-2_miR.txt")
	TGAACT<-read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/other_analyses/u_analysis/miR_analysis_3/TGAACT-2_miR.txt")
	CTAGTC<-read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/other_analyses/u_analysis/miR_analysis_3/CTAGTC-2_miR.txt")
	ATCGAA<-read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/other_analyses/u_analysis/miR_analysis_3/ATCGAA-2_miR.txt")
	GCAATT<-read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/other_analyses/u_analysis/miR_analysis_3/GCAATT-2_miR.txt")
	CAGCGT<-read.table("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/other_analyses/u_analysis/miR_analysis_3/CAGCGT-2_miR.txt")}

	
	TAGTGC<-TAGTGC[,c(1,6)]
	GCTACA<-GCTACA[,c(1,6)]
	AATCCG<-AATCCG[,c(1,6)]
	ACTGGT<-ACTGGT[,c(1,6)]
	TGTCAC<-TGTCAC[,c(1,6)]
	CGGTTA<-CGGTTA[,c(1,6)]
	GTCTAG<-GTCTAG[,c(1,6)]
	TGAACT<-TGAACT[,c(1,6)]
	CTAGTC<-CTAGTC[,c(1,6)]
	ATCGAA<-ATCGAA[,c(1,6)]
	GCAATT<-GCAATT[,c(1,6)]
	CAGCGT<-CAGCGT[,c(1,6)]
	
	colnames(TAGTGC)<-c("accession","u40m_m") #RPM
	colnames(GCTACA)<-c("accession","u1hr_m") #RPM
	colnames(AATCCG)<-c("accession","u2hr_m") #RPM
	colnames(ACTGGT)<-c("accession","u4hr_m") #RPM
	colnames(TGTCAC)<-c("accession","u8hr_m") #RPM
	colnames(CGGTTA)<-c("accession","uSS__m") #RPM
	colnames(GTCTAG)<-c("accession","u40m_p") #RPM
	colnames(TGAACT)<-c("accession","u1hr_p") #RPM
	colnames(CTAGTC)<-c("accession","u2hr_p") #RPM
	colnames(ATCGAA)<-c("accession","u4hr_p") #RPM
	colnames(GCAATT)<-c("accession","u8hr_p") #RPM
	colnames(CAGCGT)<-c("accession","uSS__p") #RPM
	
	m1<-join_all(list(
		TAGTGC,
		GCTACA,
		AATCCG,
		ACTGGT,
		TGTCAC,
		CGGTTA,
		GTCTAG,
		TGAACT,
		CTAGTC,
		ATCGAA,
		GCAATT,
		CAGCGT),by="accession",type="full")
	m1 <- m1[complete.cases(m1),]
	print(nrow(m1))
	print(head(m1))
	m1fc <- data.frame(
		"accession" = m1$accession,
		"u40m" =log2(m1$u40m_p/m1$u40m_m),
		"u1h" = log2(m1$u1hr_p/m1$u1hr_m),
		"u2h" = log2(m1$u2hr_p/m1$u2hr_m),
		"u4h" = log2(m1$u4hr_p/m1$u4hr_m),
		"u8h" = log2(m1$u8hr_p/m1$u8hr_m),
		"uSs" = log2(m1$uSS__p/m1$uSS__m))
	
	m1fc <- m1fc[complete.cases(m1fc),]
	m1fc <- m1fc[is.finite(rowSums(m1fc[,-1])),]

	print(nrow(m1fc))
	write.table(m1fc,file=paste0("/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/other_analyses/u_analysis/miR_analysis_3/all_genes_",miR,".txt"),quote=FALSE,sep="\t",row.names=FALSE)
	# return()


	# nosite_tpUTR <- m1fc[which(m1fc$accession %in% miR1_NoSite_tpUTR$V1),-1]/nrow(m1fc[which(m1fc$accession %in% miR1_NoSite_tpUTR$V1),-1])
	# nosite_topsite <- m1fc[which(m1fc$accession %in% miR1_NoSite_topsite$V1),-1]/nrow(m1fc[which(m1fc$accession %in% miR1_NoSite_topsite$V1),-1])

	# site <- m1fc[which(m1fc$accession %in% miR1_TpUtr$V1),-1]/nrow(m1fc[which(m1fc$accession %in% miR1_TpUtr$V1),-1])
	# topsite <- m1fc[which(m1fc$accession %in% miR1_TopSite$V1),-1]/nrow(m1fc[which(m1fc$accession %in% miR1_TopSite$V1),-1])
	
	site <- m1[which(m1$accession %in% miR1_TpUtr$V1),-1]
	tops <- m1[which(m1$accession %in% miR1_TopSite$V1),-1]
	snos <- m1[which(m1$accession %in% miR1_NoSite_tpUTR$V1),-1]
	tpns <- m1[which(m1$accession %in% miR1_NoSite_topsite$V1),-1]

	site <- colSums(site,na.rm=TRUE)
	tops <- colSums(tops,na.rm=TRUE)
	snos <- colSums(snos,na.rm=TRUE)
	tpns <- colSums(tpns,na.rm=TRUE)
	
	print(site)
	print(tops)
	print(snos)
	print(tpns)

	site <- site[1+2*(0:11)]/site[2+2*(0:11)]
	tops <- tops[1+2*(0:11)]/tops[2+2*(0:11)]
	snos <- snos[1+2*(0:11)]/snos[2+2*(0:11)]
	tpns <- tpns[1+2*(0:11)]/tpns[2+2*(0:11)]

	site <- log2(site[7:12]/site[1:6])
	tops <- log2(tops[7:12]/tops[1:6])
	snos <- log2(snos[7:12]/snos[1:6])
	tpns <- log2(tpns[7:12]/tpns[1:6])

	m1dc <- data.frame(
		df_site = c(site-snos),
		df_topsite = c(tops-tpns))
	m1dc$tp <- c("u40m",
		"u1h",
		"u2h",
		"u4h",
		"u8h",
		"uSS_")
	
	m1dc_melt <- melt(m1dc,id="tp")
	m1dc_melt$tp <- factor(m1dc_melt$tp,levels=c(
		"u40m",
		"u1h",
		"u2h",
		"u4h",
		"u8h",
		"uSS_"))
	print(head(m1dc_melt))
	if(miR=="miR155"){
		# p1<-ggplot(m1dc_melt[which(!m1dc_melt$variable=="nosite"),],aes(x=tp,y=log2(value),fill=variable))+geom_bar(stat="identity",position="dodge")+theme_tim()+theme(axis.line=element_blank())+scale_y_continuous(limits=c(-.5,2.5),breaks=c(-0.5,0,0.5,1.0,1.5,2.0),name=NULL,expand = c(0,0))+geom_segment(aes(y=-.5,yend=2,x=-Inf,xend=-Inf))
		p1<-ggplot(m1dc_melt[which(m1dc_melt$variable=="df_site"),],aes(x=tp,y=value))+geom_bar(stat="identity",position="dodge",fill="red")+theme_tim()+theme(axis.line=element_blank())+scale_y_continuous(limits=c(-.5,2.5),breaks=c(-0.5,0,0.5,1.0,1.5,2.0),name=NULL,expand = c(0,0))+geom_segment(aes(y=-.5,yend=2,x=-Inf,xend=-Inf))+facet_wrap(~variable)
		p2<-ggplot(m1dc_melt[which(m1dc_melt$variable=="df_topsite"),],aes(x=tp,y=value))+geom_bar(stat="identity",position="dodge",fill="blue")+theme_tim()+theme(axis.line=element_blank())+scale_y_continuous(limits=c(-.5,2.5),breaks=c(-0.5,0,0.5,1.0,1.5,2.0),name=NULL,expand = c(0,0))+geom_segment(aes(y=-.5,yend=2,x=-Inf,xend=-Inf))+facet_wrap(~variable)
		ggsave(plot=p1,file="/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/figures/other/FigS14Cpanel1.pdf",width=1.5,height=2)
		ggsave(plot=p2,file="/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/figures/other/FigS14Cpanel2.pdf",width=1.5,height=2)

	}

	if(miR=="miR1"){
		# p1<-ggplot(m1dc_melt[which(!m1dc_melt$variable=="nosite" & !m1dc_melt$tp=="u8hr"),],aes(x=tp,y=log2(value),fill=variable))+geom_bar(stat="identity",position="dodge")+theme_tim()+theme(axis.line=element_blank())+scale_y_continuous(limits=c(-.5,2.5),breaks=c(-0.5,0,0.5,1.0,1.5),name=NULL,expand = c(0,0))+geom_segment(aes(y=-.5,yend=1.5,x=-Inf,xend=-Inf))
		p1<-ggplot(m1dc_melt[which(m1dc_melt$variable=="df_site" & !m1dc_melt$tp=="u8h"),],aes(x=tp,y=value))+geom_bar(stat="identity",position="dodge",fill="red")+theme_tim()+theme(axis.line=element_blank())+scale_y_continuous(limits=c(-.5,2.5),breaks=c(-0.5,0,0.5,1.0,1.5,2.0),name=NULL,expand = c(0,0))+geom_segment(aes(y=-.5,yend=2,x=-Inf,xend=-Inf))+facet_wrap(~variable)
		p2<-ggplot(m1dc_melt[which(m1dc_melt$variable=="df_topsite" & !m1dc_melt$tp=="u8h"),],aes(x=tp,y=value))+geom_bar(stat="identity",position="dodge",fill="blue")+theme_tim()+theme(axis.line=element_blank())+scale_y_continuous(limits=c(-.5,2.5),breaks=c(-0.5,0,0.5,1.0,1.5,2.0),name=NULL,expand = c(0,0))+geom_segment(aes(y=-.5,yend=2,x=-Inf,xend=-Inf))+facet_wrap(~variable)
		ggsave(plot=p1,file="/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/figures/other/FigS14Bpanel1.pdf",width=1.5,height=2)
		ggsave(plot=p2,file="/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/figures/other/FigS14Bpanel2.pdf",width=1.5,height=2)
	}
}
plot_miR("miR1")
plot_miR("miR155")