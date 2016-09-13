# Extracting reads mapping to BRCA1 from NA12878 HiSeq Exome dataset  

1. Download the reference genome file  
`$ wget http://vannberg.biology.gatech.edu/data/ahcg2016/reference_genome/hg19_refGene.txt`  
  
2. Grab the genome coordinates for BRCA1  
`$ grep -w "NM_007294" hg19_refGene.txt > NM_007294.txt`  
&nbsp;&nbsp;Note: NM_007294 is the popular transcript of BRCA1  

3. Extract chrom, chromStart, chromEnd information from NM_007294 and convert it to bed file format.  
`$ awk '{print $3, $5-1000, $6+1000}' NM_007294.txt > NM_007294_gene_coordinates.bed` 

4. Download the NA12878 HiSeq Exome dataset BAM files from the [link](http://ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/)    

5. Using samtools to subset the bam file to regions corresponding to BRCA1.  
`$ samtools view -L <bed file> -b -o < outout bam file > < input bam file >`  
&nbsp;&nbsp;&nbsp;&nbsp;Note: -b just specifies that the output needs to be a bam file.   
> Example code:
`$ samtools view -L NM_007294_gene_coordinates.bed -b -o BRCA1_1.bam project.NIST_NIST7035_H7AP8ADXX_TAAGGCGA_1_NA12878.bwa.markDuplicates.bam`  
`$ samtools view -L NM_007294_gene_coordinates.bed -b -o BRCA1_2.bam project.NIST_NIST7035_H7AP8ADXX_TAAGGCGA_2_NA12878.bwa.markDuplicates.bam`   
`$ samtools view -L NM_007294_gene_coordinates.bed -b -o BRCA1_3.bam project.NIST_NIST7086_H7AP8ADXX_TAAGGCGA_1_NA12878.bwa.markDuplicates.bam`   
`$ samtools view -L NM_007294_gene_coordinates.bed -b -o BRCA1_4.bam project.NIST_NIST7086_H7AP8ADXX_TAAGGCGA_2_NA12878.bwa.markDuplicates.bam` 
6. Merge bam files  
&nbsp;&nbsp;Note:Put all your bam file in a folder (or create a folder with symbolic links to all sam files you want merge) and then  
`$ samtools merge finalBamFile.bam *.bam`
  
7. Using bedtools to convert the bam file to a fastq file    
&nbsp;&nbsp;Note:From the BRCA1 bam file we now extract the reads aligning to the region using bedtools      
&nbsp;&nbsp;Note:If your BAM alignments are from paired-end sequence data, one can use the -fq2 option to create two distinct FASTQ output files — one for end 1 and one for end 2. When using this option, it is required that the BAM file is sorted/grouped by the read name. This keeps the resulting records in the two output FASTQ files in the same order.   
`$ samtools sort -n aln.bam aln.qsort.bam`    
`$ bedtools bamtofastq -i <bam file> -fq < fastq r1> -fq2 < fastq r2>`  

8. Picard genome dictionary  
&nbsp;&nbsp;Note: The reference file could be download from (link)[http://vannberg.biology.gatech.edu/data/chr17.fa]  
`$ java -jar picard.jar CreateSequenceDictionary R=ref.fa O=ref.dict`

9. Samtools fasta index   
`$ samtools faidx <path to ref.fa>`

10. Bowtie index  
`$ bowtie2-build <path to ref.fa> <output prefix>`

11. Run the ahcg\_pipeline.py to get BRCA1 variants call  
&nbsp;&nbsp;Note: More information about ahcg_pipeline could be found from [link](https://github.com/mnzhang2012/ahcg_pipeline/blob/master/  README.md#variant-calling-pipeline-for-genomic-data-analysis)  
&nbsp;&nbsp;Note: here is the example command.  
` $ python3 ahcg_pipeline.py -t lib/Trimmomatic-0.36/trimmomatic-0.36.jar -b lib/bowtie2-2.2.9/bowtie2 -p lib/picard.jar -g lib/GenomeAnalysisTK.jar -i lab/BRCA1_bam/merge_BRCA1_end1.fq lab/BRCA1_bam/merge_BRCA1_end2.fq -w lab/chr17 -d dbsnp/dbsnp_138.hg19.vcf -r lab/chr17.fa -a lib/Trimmomatic-0.36/adapters/NexteraPE-PE.fa -o output_merge_BRC1`  
&nbsp;&nbsp;Note: The output file could be found from vcf_files/merge_BRCA1.vcf