#! /usr/local/env python3.4
#argv[1] is the vcf file 
#argv[2] is the reads depth information file 
#argv[3] is the final vcf for the individual

import vcf
import sys

print("# This vcf file represents all the dcm-related variants in this individual. Each variant has the gene_name, exon_number, and reads_depth information")
print("CHROM" + "\t" + "POS" + "\t" + "REF" + "\t" + "ALT" + "\t" + "Gene_Name" + "\t" + "Which_exon_it_located"  + "\t" + "Depth" + "\t" + "QUAL" + "\t" + "FILTER" + "\t" + "INFO" + "\t" + "FORMAT")
	
vcff = open(sys.argv[1])

with open(sys.argv[2]) as fin:
	lines = fin.readlines()
lines.pop(0)
with open(sys.argv[3],"w") as output:
	output.write("# This vcf file is the final vcf for all the dcm-related variants in this individual. The variants listed here are the variants passed gatk-recalibrator filter and also has the read_depth number greater than 20" + "\n")
	output.write("CHROM" + "\t" + "POS" + "\t" + "REF" + "\t" + "ALT" + "\t" + "Gene_Name" + "\t" + "Which_exon_it_located"  + "\t" + "Depth" + "\t" + "QUAL" + "\t" + "FILTER" + "\t" + "INFO" + "\t" + "FORMAT" + "\n")

	for record in vcf.Reader(vcff):
		for line in lines:
			position = line.split()[1].strip()
			gene = line.split()[2].strip()
			exon = line.split()[3].strip()
			depth = line.split()[4].strip()
			if record.POS == int(position):
				if record.FILTER is None:
					filterVar = "."
				elif len(record.FILTER) > 0:
					filterVar = str(record.FILTER[0])
				else:
					filterVar = "PASS"
					if int(depth) > 20:
						output.write(str(record.CHROM) + "\t" + str(record.POS) + "\t" + str(record.REF) + "\t" + str(record.ALT[0]) + "\t" + gene + "\t" + exon + "\t" + depth + "\t" + str(record.QUAL) + "\t" + filterVar + "\t" + str(record.INFO) + "\t" + str(record.FORMAT) + "\n")
				print(str(record.CHROM) + "\t" + str(record.POS) + "\t" + str(record.REF) + "\t" + str(record.ALT[0]) + "\t" + gene + "\t" + exon + "\t" + depth + "\t" + str(record.QUAL) + "\t" + filterVar + "\t" + str(record.INFO) + "\t" + str(record.FORMAT))

