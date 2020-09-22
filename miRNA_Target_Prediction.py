import re
import time
import sys
import os
import pandas as pd
import requests
import random

def getArgvDict(argv):
		optionDict = {}
		c = 0
		for i in argv:
				if re.match(r'-\D', i):
						optionDict[i] = argv[c + 1]
				c += 1
		return optionDict

def ScrapyTarScan(url, regEx):
		time.sleep(random.uniform(2.2,4.4))
		response = requests.get(url)
		targetscanFile = response.content
		targetscanResPage = targetscanFile.decode('utf-8')
		regObj = re.compile(regEx, re.DOTALL)
		match = regObj.search(targetscanResPage)
		if match:
				scrapyRes = match.group(1)
				return scrapyRes
		else:
				return 'Nothing'

def GetTargetScan(miR, out, threshold):
		baseUrl1 = 'http://www.targetscan.org/cgi-bin/targetscan/vert_72/'
		baseUrl2 = 'http://www.targetscan.org'
		searchItem = 'targetscan.cgi?species=Human&gid=&mir_sc=&mir_c=&mir_nc=&mir_vnc=&mirg='
		regEx1 = '<A HREF="(.*?)">Download'
		regEx2 = '<A HREF="(.*?)" download'
		regEx3 = '(.*)'
		targetscanUrl = baseUrl1 + searchItem + miR
		downloadItem = ScrapyTarScan(targetscanUrl, regEx1)
		if 'Nothing' == downloadItem:
				print('Warning: This miRNA('+miR+') is not exists in TargetScan v7.2, you may need check this out mannully.')
				geneList = []
		else:
				downloadUrl = baseUrl1 + downloadItem
				textItem = ScrapyTarScan(downloadUrl, regEx2)
				textUrl = baseUrl2 + textItem
				targetRes = ScrapyTarScan(textUrl, regEx3)
				t = open(out + 'TargetScan_v7.2/TMP', 'w')
				t.write(targetRes)
				t.close()
				oriTSC = pd.read_table(out + 'TargetScan_v7.2/TMP', header = 0, sep = '\t')
				tmp1 = oriTSC[oriTSC['Representative miRNA'] == miR]
				tmp2 = tmp1[tmp1['Cumulative weighted context++ score'] <= float(threshold)]
				geneList = tmp2['Target gene'].values
				tmp2.to_csv(out + 'TargetScan_v7.2/TarScan-' + miR + '.tsv', index = False, header = True, sep = '\t')
				os.remove(outDir + 'TargetScan_v7.2/TMP')
		return geneList

def GetMiRDB(miR, RDB, out, threshold):
		oriRDB = pd.read_table(RDB, header = 0, sep = '\t')
		tmp1 = oriRDB[oriRDB['miRNA'] == miR]
		if len(tmp1) <= 1:
				print('Warning: This miRNA('+miR+') is not exists in miRDB v6.0, you may need check this out mannully.')
				geneList = []
		else:
				tmp2 = tmp1[tmp1['score'] >= int(threshold)]
				geneList = tmp2['Gene'].values
				tmp2.to_csv(out + 'miRDB_v6.0/miRDB_NM-' + miR + '.tsv', index = False, header = True, sep = '\t')
		return geneList

def GetFinalTarget(miR, tsc, rdb, out):
		clipdf = pd.read_table('./data/GSE111432_AGO2CLIP_gene.txt', header = None, names=['SYMBOL'], sep = '\t')
		clipSymbol = clipdf['SYMBOL'].values
		with open('./data/miRDB_v6.0.hsaOnly.RefSeq2SYMBOL.txt', 'r') as f:
				cliptxt = f.read()
		rdb_symbol = []
		for i in rdb:
				match = re.search(i+r'\t(\S+)\n', cliptxt)
				if match:
						rdb_symbol.append(match.group(1))
		TSC_b_RDB = set(tsc).union(set(rdb_symbol))
		TSC_b_RDB_j_CLIP = TSC_b_RDB.intersection(set(clipSymbol))
		clip = [s for s in TSC_b_RDB_j_CLIP if s != '']
		t = open(out+'final_targets/Targets-' + miR + '.tsv', 'w')
		t.writelines("\n".join(TSC_b_RDB_j_CLIP))
		t.close()




if __name__ == '__main__':
		helpDoc = '''
-H/-h	Show this help-doc.

-miRNA	The miRNAs list file , one miRNA per line. Please use the
	complete miRNA name, e.g. hsa-miR-218-5p.

[-tarScanScore]	Optional. The threshold of cumulative weighted 
	context++ score in TargetScan-V7.2 to filter results. The
	predicted target with the score great than this threshold
	will be filtered out. The value range is [-1, 0] (-1 means
	best prediction). The defult value is -0.2.

[-miRDBScore]	Optional. The threshold for prediction confidence
	in miRDB-V6.0 to filter results. The predicted target with
	the score lower than this threshold will be filtered out. 
	The value range is [50, 100] (100 means best prediction).
	The defult value is 60.

-outputDir	A folder will be used to store results.

The script will output the following folder/files into -outputDir:

	TargetScan_v7.2/TarScan-miRNA.tsv, TargetScan results.
	miRDB_v6.0/miRDB_NM-miRNA.tsv, miRDB results.
	final_targets/Targets-miRNA.tsv, the union of TargetScan and
	miRDB results, and filterd by AGO2 CLIP-seq data (osteoblast).
		'''
		argv = sys.argv
		if re.search('-H|-h', str(argv)):
				print(helpDoc)
				sys.exit()
		try:
				envDict = getArgvDict(argv)
				miRFile = envDict['-miRNA']
				outDir = envDict['-outputDir']
		except KeyError as e:
				print('ERROR: Incomplete or invalid arg '+str(e)+' ! Please check your input, or use -H/-h flag to get help.')
				sys.exit()
		else:
				tarScanScore = (envDict['-tarScanScore'] if envDict.get('-tarScanScore') else -0.2)
				miRDBScore = (envDict['-miRDBScore'] if envDict.get('-miRDBScore') else 60)
				if not os.path.exists(outDir):
						os.makedirs(outDir)
				outDir = (outDir if outDir[-1] == '/' else outDir + '/')
				if not os.path.exists(outDir+'TargetScan_v7.2'):
							os.makedirs(outDir+'TargetScan_v7.2')
				if not os.path.exists(outDir+'miRDB_v6.0'):
							os.makedirs(outDir+'miRDB_v6.0')
				if not os.path.exists(outDir+'final_targets'):
							os.makedirs(outDir+'final_targets')
				miRDBFile = './data/miRDB_v6.0_hsaOnly.txt'
		with open(miRFile, 'r') as f:
				miRList = f.readlines()
		for miR in miRList:
				miR = miR.strip()
				print('Predicting the targets of ' + str(miR) + '...')
				tarListTSC = GetTargetScan(miR, outDir, tarScanScore)
				tarListRDB = GetMiRDB(miR, miRDBFile, outDir, miRDBScore)
				GetFinalTarget(miR, tarListTSC, tarListRDB, outDir)
		localtime = time.asctime(time.localtime(time.time()))
		print('Done. Finish time: '+ localtime + '\n')


