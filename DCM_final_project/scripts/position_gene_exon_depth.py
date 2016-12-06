import sys
#argv[1] is the file which contains two columns. The first column contains name of gene. The 2nd column contains the NM_. The file has header.
#argv[2] is human reference genome file, which you could extract the genome coordinates by given the NM_ number.
#argv[3] is the final depth file for each patient

with open(sys.argv[1]) as fin:
	match = []
	for line_fin in fin:
		with open(sys.argv[2]) as fin_ref:
			for line_ref in fin_ref:
				if line_fin.split()[1].rstrip() in line_ref:
					match.append(line_ref + "\t" + line_fin.split()[0].rstrip())				

with open("gene_list_BED.txt","w") as output_bed:
	for line in match:
		chrom_start = line.split()[9].split(",")
		chrom_end = line.split()[10].split(",")
		for x in range(0,len(chrom_start)-1):
			result = line.split()[2]+"\t"+chrom_start[x]+"\t"+chrom_end[x]+"\t"+line.split()[-1]
			output_bed.write(result + "\n")

with open ("gene_list_BED.txt") as input_bed:
	with open ("gene_list_BED_add_20.txt","w") as output_bed:
		for line in input_bed:
			line20 = str(line.split()[0] + "\t" + str(int(line.split()[1]) - 20) + "\t" + str(int(line.split()[2]) + 20) + "\t" + line.split()[-1])
			output_bed.write(line20 + "\n")


with open ("gene_list_BED_add_20.txt") as fin:
	lines = fin.readlines()

with open("gene_list_BED_add_20_exon.txt","w") as output_bed:
	first = ""
	for line in lines:
		gene = line.split()[3].strip()
		if gene == first:
			newline = str(line.strip()) + "\t" + str(count)
			output_bed.write(newline + "\n")
		else:
			count = 1
			newline = str(line.strip()) + "\t" + str(count)
			output_bed.write(newline + "\n")
		first = gene
		count = count + 1

total = 0
with open ("gene_list_BED_add_20_exon.txt") as f:
	with open("gene_list_BED_add_20_exon_region_to_position.txt","w") as output_bed:
		for line in f:
			array = line.split()
			chrom = array[0]
			start = int(array[1])
			end = int(array[2])
			gene = array[3]
			exon = array[4]
			while (start < end):
				output_bed.write(chrom + "\t" + str(start) + "\t" + gene  + "\t" + exon + "\n")
				start = start + 1

with open("gene_list_BED_add_20_exon_region_to_position.txt") as fin1:
	exons = fin1.readlines()

with open(sys.argv[3]) as fin2:
	depths = fin2.readlines()

def header():
	print("CRHOM"+"\t"+"POSITION"+"\t"+"GENE"+"\t"+"Which_exon_it_located"+"\t"+"Read_Depth")

header()

for exon in exons:
	for depth in depths:
		if exon.split()[1] == depth.split()[1]:
			print(exon.strip() + "\t" + depth.split()[-1].strip())