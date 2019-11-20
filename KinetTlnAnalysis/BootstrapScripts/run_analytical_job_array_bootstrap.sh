#!/bin/bash

BootstrapNo=1

while [ $BootstrapNo -le 2 ]
do
	file=model_input_files/bootstrap/miR-155_minus_sample_mean_tails_50tags_with_expr_background_subtracted_v6_st_BOOTSTRAP_${BootstrapNo}.txt
	N=$(($(wc -l < $file)-1)) #remove the last line because of the header
	#N=2
	time=$((BootstrapNo-1+22))
	bsub -b $time:00 bsub -o /lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/processed_files/bootstrap/bsub_log.txt -J "V${BootstrapNo}[1-$N]%280" "Rscript /lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/KinetTlnAnalysis/Model/AnalyticalTailsPulseStUnlinkedBootstrap.R $file \$LSB_JOBINDEX v$BootstrapNo $BootstrapNo"
	((BootstrapNo++))
done
while [ $BootstrapNo -le 10 ]
do
	file=model_input_files/bootstrap/miR-155_minus_sample_mean_tails_50tags_with_expr_background_subtracted_v6_st_BOOTSTRAP_${BootstrapNo}.txt
	N=$(($(wc -l < $file)-1)) #remove the last line because of the header
	#N=2
	time=$((BootstrapNo-1-2))
	bsub -b $time:00 bsub -o /lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/processed_files/bootstrap/bsub_log.txt -J "V${BootstrapNo}[1-$N]%280" "Rscript /lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/KinetTlnAnalysis/Model/AnalyticalTailsPulseStUnlinkedBootstrap.R $file \$LSB_JOBINDEX v$BootstrapNo $BootstrapNo"
	((BootstrapNo++))
done
# while [ $BootstrapNo -le 10 ]
# do
# #needs to bsub with dependencies
# 	file=model_input_files/bootstrap/miR-155_minus_sample_mean_tails_50tags_with_expr_background_subtracted_v5_st_BOOTSTRAP_${BootstrapNo}.txt
# 	# N=$(($(wc -l < $file)-1)) #remove the last line because of the header
# 	N=1
# 	bsub -o /lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/processed_files/bootstrap/bsub_log.txt -J "V${BootstrapNo}[1-$N]%280" -w "ended(\"V$((BootstrapNo-1))\")" "echo /lab/solexa_bartel/teisen/RNAseq/Scripts/models/single_tail_fitting/analytical_tails_V75_pulse_st.R $file \$LSB_JOBINDEX v$BootstrapNo $BootstrapNo"
# 	((BootstrapNo++))
# done
