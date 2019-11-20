#!/bin/bash

# COUNTER=1
# while [  $COUNTER -lt 11 ]; do
# 	bsub -q 18 -o /lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/rate_constant_measurements/global_model_fits/bootstrap/run_output_HYBRID.txt -J "HybridBootstrap$COUNTER[1-30]" "Rscript /lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/KinetTlnAnalysis/Model/UnLinkV3/AnalyticalTailsPulseStUnlinkedGlobalHybridBootstrapV4_20.R \$LSB_JOBINDEX $COUNTER"
# 	let COUNTER=COUNTER+1 
# done

# while [  $COUNTER -lt 11 ]; do
# 	bsub -b 28:20:00 -o /lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/rate_constant_measurements/global_model_fits/bootstrap/run_output.txt -J "V$COUNTER[1-30]" "Rscript /lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/KinetTlnAnalysis/Model/AnalyticalTailsPulseStUnlinkedGlobalOneDcpBootstrap.R \$LSB_JOBINDEX $COUNTER"
# 	let COUNTER=COUNTER+1 
# done

# bsub -q 18 -R "rusage[mem=4096]" -J "VGlobalHybrid$COUNTER[1-30]" "Rscript /lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/KinetTlnAnalysis/Model/UnLinkV3/AnalyticalTailsPulseStUnlinkedGlobalHybridV4_20.R \$LSB_JOBINDEX"

# bsub -q 18 -R "rusage[mem=4096]" -J "VGlobalHybridmiR1$COUNTER[1-30]" "Rscript /lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/KinetTlnAnalysis/Model/UnLinkV3/AnalyticalTailsPulseStUnlinkedGlobalHybridmiR1V4_20.R \$LSB_JOBINDEX"

#global two dea models
# bsub -q 18 -R "rusage[mem=4096]" -J "VGlobalHybrid$COUNTER[1-30]" "Rscript /lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/KinetTlnAnalysis/Model/UnLinkV3/AnalyticalTailsPulseStUnlinkedGlobalTwoDeaHybrid110V4.R \$LSB_JOBINDEX"

# bsub -q 18 -R "rusage[mem=4096]" -J "VGlobalHybrid$COUNTER[1-30]" "Rscript /lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/KinetTlnAnalysis/Model/UnLinkV3/AnalyticalTailsPulseStUnlinkedGlobalTwoDeaHybrid150V4.R \$LSB_JOBINDEX"