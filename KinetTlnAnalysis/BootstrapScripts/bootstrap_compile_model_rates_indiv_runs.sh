#!/usr/bin/bash

DatasetNo=1
while [ $DatasetNo -le 10 ]
do
	cat optim_runs/v${DatasetNo}/*.txt > rate_constant_measurements/bootstrap/miR-155_minus_miR-155_minus_samples_UNLINKV3_V81_H5_Run4_BOOTSTRAP_${DatasetNo}.txt
	((DatasetNo++))
done
