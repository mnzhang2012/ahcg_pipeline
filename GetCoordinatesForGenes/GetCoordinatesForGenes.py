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
	for line in infile:
		print(line.split()[2]+"\t"+line.split()[4]+"\t"+line.split()[5]+"\t"+line.split()[1]+"\n")
		outfile.write(line.split()[2]+"\t"+line.split()[4]+"\t"+line.split()[5]+"\t"+line.split()[1]+"\n")
infile.close()
outfile.close()
os.remove("otogenetics_BRC_biomarker.txt")