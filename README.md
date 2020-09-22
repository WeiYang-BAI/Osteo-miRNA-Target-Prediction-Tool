# Osteo-miRNA-Target-Prediction-Tool

This is an easily-to-use tool for users to predict human osteoblast-specific target genes of human miRNA of interest. The TargetScan v7.2 and miRDB v6.0 are jointly used to form a putative targets pool with the given thresholds, the AGO2 CLIP-seq data is used to fliter false postive targets. User could just input a list of miRNA in local.

Use the follow command to download:

	$ git clone https://github.com/WeiYang-BAI/Osteo-miRNA-Target-Prediction-Tool.git

The miRDB data used in the tool was gziped, please unzip it before use:

	$ cd data
	$ gunzip miRDB_v6.0_hsaOnly.txt.gz

Run the follow command to see the help message:

	$ cd ..
	$ python ./miRNA_Target_Prediction.py -H

Arguments:

	-H/-h	Show this help-doc.

	-miRNA	The miRNAs list file , one miRNA per line. Please use the complete miRNA name, e.g. hsa-miR-218-5p.

	[-tarScanScore]	Optional. The threshold of cumulative weighted context++ score in TargetScan v7.2 to filter 
			results. The predicted target with the score great than this threshold will be filtered out.
			The value range is [-1, 0] (-1 means best prediction). The defult value is -0.2.

	[-miRDBScore]	Optional. The threshold for prediction confidence in miRDB-V6.0 to filter results. The predicted
			target with the score lower than this threshold will be filtered out. The value range is [50, 100]
			(100 means best prediction). The defult value is 60.

	-outputDir	A folder will be used to store results.


An example is given, run it in the command line for details:

	$ sh example.sh


Output folder/files into -outputDir:

	TargetScan_v7.2/TarScan-miRNA.tsv, prediction of TargetScan results.
	miRDB_v6.0/miRDB_NM-miRNA.tsv, prediction of miRDB results.
	final_targets/Targets-miRNA.tsv, the union of TargetScan and miRDB results, and filterd by AGO2 CLIP-seq data.
