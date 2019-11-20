#!/usr/bin/bash

rm temp.txt
BootstrapNo=1
while [ $BootstrapNo -le 10 ]
do
	DatasetNo=1
	while [ $DatasetNo -le 34 ]
	do
		cut -f 5-8 rate_constant_measurements/global_model_fits/bootstrapHYBRID/miR-155_minus_rand_g_rates_plogis_unlinked_HYBRID20190801_${DatasetNo}_BOOTSTRAP_${BootstrapNo}.txt | head -2 | tail -n +2 | awk -v dat="$DatasetNo" -v bs="$BootstrapNo" 'BEGIN{FS=OFS="\t"}{print bs,dat,$0}' >> temp.txt
		((DatasetNo++))
	done
	((BootstrapNo++))
done

Rscript KinetTlnAnalysis/BootstrapScripts/compile_global_rates.R temp.txt /lab/solexa_bartel/teisen/Tail-seq/miR-155_final_analyses/rate_constant_measurements/global_model_fits/bootstrapHYBRID/miR-155_minus_rand_g_rates_plogis_unlinked_HYBRID20190831_compiled_summarize.txt

rm temp.txt