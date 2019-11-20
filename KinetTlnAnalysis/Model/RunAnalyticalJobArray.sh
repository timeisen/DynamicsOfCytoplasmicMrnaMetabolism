#!/bin/bash

N=$(($(wc -l < $1)-1)) #remove the last line because of the header
bsub -q 18 -o /lab/solexa_bartel/teisen/bsub_output_6.txt -J "V75[1-$N]%280" "Rscript KinetTlnAnalysis/Model/UnLinkV3/AnalyticalTailsPulseStUnlinkedHYBRIDV5RandI.R $1 \$LSB_JOBINDEX $2"

# bsub -q 18 -b 15:00 bsub -q 18 -o /lab/solexa_bartel/teisen/bsub_output_6.txt -J "V75[1-2778]%280" "Rscript KinetTlnAnalysis/Model/UnLinkV3/AnalyticalTailsPulseStUnlinkedTwoDea110V5.R model_input_files/miR-155_minus_sample_mean_tails_50tags_with_expr_background_subtracted_v7_st_last20rm_HYBRID20190731.txt \$LSB_JOBINDEX v22"

# bsub -q 18 -b 17:00 bsub -q 18 -o /lab/solexa_bartel/teisen/bsub_output_6.txt -J "V75[1-2778]%280" "Rscript KinetTlnAnalysis/Model/UnLinkV3/AnalyticalTailsPulseStUnlinkedTwoDea150V5.R model_input_files/miR-155_minus_sample_mean_tails_50tags_with_expr_background_subtracted_v7_st_last20rm_HYBRID20190731.txt \$LSB_JOBINDEX v30"

# 13.43797
# 266.8899
# 9.499726


# colnames(mm) <- c("accession","stm","am","km","bm","rm")
# colnames(mp) <- c("accession","kp","bp","rp")
# m <- merge(mm,mp,by="accession")
# m2 <- m
# m2$kl2fc <- log2(m2$kp/m2$km)
# m2$bl2fc <- log2(m2$bp/m2$bm)
# m2 <- m2[,c(1,10,11)]