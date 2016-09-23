import re
import os
with open("otogenetics_BRC_list.txt") as inf_bed:
	next(inf_bed)
	outf = open("otogenetics_BRC_biomarker.txt","w")
	for line_bed in inf_bed:
		with open("hg19_refGene.txt") as inf_ref:
			for line_ref in inf_ref:
				if re.search(line_bed.split()[1].rstrip(),line_ref):
					outf.write(line_ref)
outf.close()
inf_ref.close()
inf_bed.close()

with open("otogenetics_BRC_biomarker.txt","r") as infile, open("otogenetics_BRC_biomarker_BED.txt","w") as outfile:
	outfile.write("chr"+"\t"+"chrom_start"+"\t"+"chrom_end"+"\t"+"strand"+"\t"+"name"+"\n")
	for line in infile:
		chrom_start = line.split()[9].split(",")
		del chrom_start[-1]
		chrom_end = line.split()[10].split(",")
		del chrom_end[-1]
		for x in range(0,len(chrom_start)):
				outfile.write(line.split()[2]+"\t"+chrom_start[x]+"\t"+chrom_end[x]+"\t"+line.split()[3]+"\t"+line.split()[1]+"\n")
infile.close()
outfile.close()
os.remove("otogenetics_BRC_biomarker.txt")