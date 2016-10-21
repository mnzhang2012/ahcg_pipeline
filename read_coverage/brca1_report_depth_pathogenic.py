# argv[1] is the pathogenic annotation of variants call
# argv[2] is the variants read depth file
import sys
pathogenics = []
depths = []
with open(sys.argv[1]) as fin_pathogenic:
	with open(sys.argv[2]) as fin_depth:
		with open("brca_depth_pathogenic_report.txt","w") as fout:
			pathogenics = fin_pathogenic.readlines()
			depths = fin_depth.readlines()

			fout.write("CHROM" + "\t" + "POSITION" + "\t" + "REF" + "\t" + "ALT" + "\t" + "PATHOGENIC" + "\n")
			for pathogenic in pathogenics:
				find = 0
				for depth in depths:
					if pathogenic.split(",")[1][5:].strip() == depth.split()[1].strip():						
						find = 1
						depth_number = depth.split()[2].strip()
			
				line = pathogenic.split(",")[0][-5:] + "\t" + pathogenic.split(",")[1][5:] + "\t" + pathogenic.split(",")[2][5:] + "\t" + pathogenic.split(",")[3][6] + "\t" + pathogenic.split(")")[1].strip()
				if find == 1:
					fout.write(line.strip() + "\t" + depth_number + "\n")
				else:
					fout.write(line.strip() + "\t" + "NA" + "\n")