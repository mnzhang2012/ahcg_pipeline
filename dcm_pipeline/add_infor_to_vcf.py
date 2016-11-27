#! /usr/local/env python3.4
#argv[1] is the vcf file 

import vcf
import sys

vcff  =  open(sys.argv[1])
with open(sys.argv[2]) as fin:
	lines = fin.readlines()

print("CHROM" + "\t" + "POS" + "\t" + "REF" + "\t" + "ALT" + "\t" + "QUAL" + "\t" + "FILTER" + "\t" + "INFO" + "\t" + "FORMAT" + "\t" + "DEPTH")
	
for record in vcf.Reader(vcff):
	for line in lines:
		position = line.split()[1].strip()
		depth = line.split()[2].strip()
		if record.POS == int(position):
			print(str(record.CHROM) + "\t" + str(record.POS) + "\t" + str(record.REF) + "\t" + str(record.ALT) + "\t" + str(record.QUAL) + "\t" + str(record.FILTER) + "\t" + str(record.INFO) + "\t" + str(record.FORMAT) + "\t" + depth)