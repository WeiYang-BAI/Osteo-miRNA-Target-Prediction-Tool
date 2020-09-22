#!/bin/sh

python \
	./miRNA_Target_Prediction.py \
	-miRNA ./example/example-miR.txt \
	-tarScanScore -0.2 \
	-miRDBScore 60 \
	-outputDir ./example/prediction_results
