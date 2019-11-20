#!/bin/bash

########################################################
# run_analytical_job_array_bootstrap_HYBRID.sh
# This code runs the bootstrap models for individual gene fitting.
# Make sure that the Rscript is the bootstrap version of the individual tail fitting
# It must be able to take global input params
# TJE 2019 05 15, last update.
# See github for version history. 
# 
########################################################
BootstrapNo=1

# while [ $BootstrapNo -le 2 ]
# do
# 	file=model_input_files/bootstrap/miR-155_minus_sample_mean_tails_50tags_with_expr_background_subtracted_v7_st_BOOTSTRAP_HYBRID20190801_Rescale_${BootstrapNo}.txt
# 	N=$(($(wc -l < $file)-1)) #remove the last line because of the header
# 	#N=2
# 	time=$((BootstrapNo-1+22))
# 	bsub -q 18 -b $time:00 bsub -q 18 -o /lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/processed_files/bootstrap/bsub_log.txt -J "V${BootstrapNo}[1-$N]%280" "Rscript /lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/KinetTlnAnalysis/Model/UnLinkV3/AnalyticalTailsPulseStUnlinkedBootstrapHYBRIDV5.R $file \$LSB_JOBINDEX v$BootstrapNo $BootstrapNo"
# 	((BootstrapNo++))
# done
while [ $BootstrapNo -le 10 ]
do
	file=model_input_files/bootstrap/miR-155_minus_sample_mean_tails_50tags_with_expr_background_subtracted_v7_st_BOOTSTRAP_HYBRID20190801_Rescale_${BootstrapNo}.txt
	N=$(($(wc -l < $file)-1)) #remove the last line because of the header
	#N=2
	time=$((BootstrapNo + 13))
	bsub -q 18 -b $time:00 bsub -q 18 -o /lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/processed_files/bootstrap/bsub_log.txt -J "V${BootstrapNo}[1-$N]%280" "Rscript /lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/KinetTlnAnalysis/Model/UnLinkV3/AnalyticalTailsPulseStUnlinkedBootstrapHYBRIDV5.R $file \$LSB_JOBINDEX v$BootstrapNo $BootstrapNo"
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

# 1-243,252-493,502-743,752-993,1002-1243,1252-1493