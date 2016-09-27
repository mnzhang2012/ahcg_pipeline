with open ("BRC_OC_gene_list_BED.txt") as input_bed:
	output_bed = open("BRC_OC_gene_list_BED_add_20.txt","w")
	next(input_bed)
	for line in input_bed:
		line20 = str(line.split()[0] + "\t" + str(int(line.split()[1]) - 20) + "\t" + str(int(line.split()[2]) + 20) + "\t" + line.split()[3] + "\t" + line.split()[4] + "\n")
		output_bed.write(line20)
input_bed.close()
output_bed.close() 