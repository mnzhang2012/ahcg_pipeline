import sys
#argv[1] is the BED file with the following format "chr	chrom_start	chrom_end	strand	name". The file has header.
with open (sys.argv[1]) as fin:
	for line in fin:
		if line.startswith("chr"):
			print(line.replace("chr",""), end="")
		else:
			print(line, end="")