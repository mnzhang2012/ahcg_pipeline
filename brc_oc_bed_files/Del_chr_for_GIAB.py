import sys
with open (sys.argv[1]) as input_bed:
	for line in input_bed:
		if line.startswith("chr"):
			print(line.replace("chr",""), end="")
		else:
			print(line, end="")
input_bed.close()
 
