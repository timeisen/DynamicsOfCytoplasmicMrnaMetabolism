#!/usr/bin/bash

BootstrapNo=1
while [ $BootstrapNo -le 10 ]
do
	DatasetNo=1
	while [ $DatasetNo -le 34 ]
	do
		cut -f 5-8 rate_constant_measurements/global_model_fits/bootstrapHYBRID/miR-155_minus_rand_g_rates_plogis_unlinked_HYBRID_${DatasetNo}_BOOTSTRAP_${BootstrapNo}.txt | head -2 | tail -n +2 | awk -v dat="$DatasetNo" -v bs="$BootstrapNo" 'BEGIN{FS=OFS="\t"}{print bs,dat,$0}'
		((DatasetNo++))
	done
	((BootstrapNo++))
done

# colnames(m) <- c("DatasetNo","bootstrap","si","loc","sc","resid")
# mtemp <- NULL
# mcompiled <- NULL
# for(i in 1:max(m$DatasetNo)){
# 	mtemp = m[m$DatasetNo == i,]
# 	mcompiled = rbind(mcompiled,apply(mtemp,2,median))
# }
# mcompiled <- data.frame(mcompiled)
# mcompiled$bootstrap <- NULL

# write.table(mcompiled,"/lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/rate_constant_measurements/global_model_fits/bootstrapHYBRID/miR-155_unlinked_rates_compiled_median_HYBRID_20180925.txt",row.names=FALSE,sep="\t",quote=FALSE)