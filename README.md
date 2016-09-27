# Variant calling pipeline for genomic data analysis

###1. Software Requirements  
1.1. Here is a list of the programs that is needed for running this pipeline  
   &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Download the programs using wget command

  1. [Python3 - version 3.4.1](https://www.python.org/download/releases/3.4.1/)
  2. [Trimmomatic - version 0.36](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip)
  3. [Bowtie2 - version 2.2.9](https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.9/)
  4. [Picard tools - version 2.6.0](https://github.com/broadinstitute/picard/releases/download/2.6.0/picard.jar)
  5. [GATK - version 3.4](https://software.broadinstitute.org/gatk/download/)  

1.2. Reference genome and Test data  
1.2.1. Reference genome
  Reference genomes can be downloaded from [Illumina iGenomes](http://support.illumina.com/sequencing/sequencing_software/igenome.html)  
 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;This pipeline is using [UCSC-hg19 reference genome]    
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`$ wget www.prism.gatech.edu/~sravishankar9/resources.tar.gz`

1.2.2. Test data  
  Use the following protocol to download and prepare test dataset from NIST sample NA12878  

 >`$ wget ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R1_001.fastq.gz`  
  `$ wget ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R2_001.fastq.gz`  
  `$ gunzip NIST7035_TAAGGCGA_L001_R1_001.fastq.gz`  
 ` $ gunzip NIST7035_TAAGGCGA_L001_R2_001.fastq.gz`  
  `$ head -1000 NIST7035_TAAGGCGA_L001_R1_001.fastq > test_r1.fastq`  
 ` $ head -1000 NIST7035_TAAGGCGA_L001_R2_001.fastq > test_r2.fastq`

###2. Instructions for running pipeline  
2.1. Build bowite2 index for reference genome  
 >`$ bowtie2-build <reference_in> <bt2_base>`  
    The prebuilt index for reference genome hg19 could be download from ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/hg19.zip  

2.2. Build fasta index for reference genome  
   `$ samtools faidx ref.fasta`

2.3. Build genome dict file for reference genome  
   `$ java -jar picard.jar CreateSequenceDictonary R=ref.fasta O=ref.dict`

2.4. Reorder the VCF file  
   `$ java -jar picard.jar SortVcf   I=original.vcf O=sorted.vcf SEQUENCE_DICTIONARY=ref.dict`

2.5. run the ahcg_pipeline.py script  
   &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Exmaple command  
   > `$ python3 ahcg_pipeline.py -t lib/Trimmomatic-0.36/trimmomatic-0.36.jar -b lib/bowtie2-2.2.9/bowtie2 -p lib/picard.jar -g lib/GenomeAnalysisTK.jar -i sequence/test_r1.fastq sequence/test_r2.fastq -w hg19 -d dbsnp/dbsnp_138.hg19.vcf -r genome/hg19.fa -a lib/Trimmomatic-0.36/adapters/NexteraPE-PE.fa -o output`   
   &nbsp;&nbsp;Notes:   
     &nbsp;&nbsp;1. -i input files option  
       &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; -i pair-end-read1 pair-end-read2  
     &nbsp;&nbsp;2. -w bowtie2 index files option  
       &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; -w hg19  
    &nbsp; 3. For the error of running out of memory for java application, you could change the give more memory for the virtual box and also make sure java's maximum memory usage has been increased according to how much you assigned for your virtual box using the arguments -Xmx<memory>.


###3. Help

To access help use the following command:
python3 ahcg_pipeline.py -h



# Readme file for change the remote url for git repositories on Mac terminal
1. Change the current working directory to local git repository folder
2. List exiting remotes for the repository  
   `$ git remote -v`
3. Use the git remote set-url command to change the remote repository's URL  
   `$ git remote set-url origin https://github.com/mnzhang2012/ahcg_pipeline.git`
4. Verify the remote URL has changed  
   `$ git remote -v`


# Readme file for updating an existing file to a GitHub repository using the command line
1. Copy the file that needed to be updated/uploaded to GitHub into the local directory that was created when you cloned the repository
2. Change the current working directory to local repository
3. Stage the file for commit to your local repository using the following command  
>`$ git add .`  
>This command adds / stages all of the files in the current directory. This is for convenience, and can still be used if you have certain files you don't want to add by using a .gitignore  

4. Commit the file that you've staged in local repository  
>`$ git commit -m "Notes related with changes"`  
>This command stores the current contents of the index in a new commit along with a log message from the user describing the changes  

5. Push the chagnes in local repository to Github remote repository  
>`$ git push remote_name branch_name`  
For example, use $ git push origin master   
Notes: origin is an alias on your system for a particular remote repository. It's not actually a property of that repository.    
Notes: A branch in Git is simply a lightweight movable pointer to one of these commits. The default branch name in Git is master. As you initially make commits, you're given a master branch that points to the last commit you made. Every time you commit, it moves forward automatically.
  
  
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
>`$ python3 ahcg_pipeline.py -t lib/Trimmomatic-0.36/trimmomatic-0.36.jar -b lib/bowtie2-2.2.9/bowtie2 -p lib/picard.jar -g lib/GenomeAnalysisTK.jar -i lab/BRCA1_bam/merge_BRCA1_end1.fq lab/BRCA1_bam/merge_BRCA1_end2.fq -w lab/chr17 -d dbsnp/dbsnp_138.hg19.vcf -r lab/chr17.fa -a lib/Trimmomatic-0.36/adapters/NexteraPE-PE.fa -o output_merge_BRC1`  
&nbsp;&nbsp;Note: The output file could be found from vcf_files/merge_BRCA1.vcf


# Compare the vcf file(running from ahcg_pipeline) with golden standard vcf file.   
1. Download the BED file, which contains the Nextera Rapid Capture Expanded Exome targeted regions.   
`$ wget http://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/samplepreps_nextera/nexterarapidcapture/nexterarapidcapture_expandedexome_targetedregions.bed`  
  
2. Download the vcf file running from ahcg_pipeline.  
`$ wget http://vannberg.biology.gatech.edu/data/variants.vcf`  
  
3. Download the golden standard vcf file.  
`$ wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/analysis/10XGenomics_calls_08142015/NA12878_phased_variants.vcf.gz`

4. Use bedtools to cut the Nextera Rapid Capture Expanded Exome targeted regions from golden standard vcf file.  
`$ bedtools intersect -wa -a <vcf file> -b <bed file>`  
  
5. Use bedtools to find intersect calls between vcf(from pipeline) and golden standard vcf.  
`$ bedtools intersect -header -a <vcf_golden_standard> -b <vcf_pipeline>`   
  
6. Use bedtools subtract to find calls that are unique for gold_standard vcf.  
`$ bedtools substract -a <vcf_golden_standard> -b <vcf_intersect>`
  
7. Use bedtools subtract to find calls that are unique for pipeline vcf.  
`$ bedtools substract -a <vcf_pipeline> -b <vcf_intersect>`
  
8. Use bedtools jaccard to find the number of TP, FP, FN variant calls.
`$ bedtools jaccard -a <vcf_intersect> -b <vcf_pipeline>`
`$ bedtools jaccard -a <vcf_intersect> -b <vcf_golden_standard>`  
  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Note: The output file could be found from vcf_files/golden_standard_variants.vcf



# Extracting chromosome coordinates given NCBI_Accession number.    

1. Download the reference genome file  
`$ wget http://vannberg.biology.gatech.edu/data/ahcg2016/reference_genome/hg19_refGene.txt`   

2. Download the list of breast cancer biomarkers.  
`$ wget http://www.otogenetics.com/forms/Breast_Cancer_gene_list.pdf`   
&nbsp;&nbsp;&nbsp;&nbsp;Note: Extract the first two columns and exported as a txt file.    
&nbsp;&nbsp;&nbsp;&nbsp;Note: The first column is the gene symbols, the second column is the NCBI reference.    

3. Run the python pipeline  
&nbsp;&nbsp;&nbsp;&nbsp;Note: The python pipeline is in the folder GetCoordinatesForGenes.  
`$ python GetCoordinatesForGenes.py`  
&nbsp;&nbsp;&nbsp;&nbsp;Note: The example output is in the folder GetCoordinatesForGenes.  

# Compare the ahcg\_pipeline vcf file with GIAB vcf file for genes from color genomics and otogenetics, breast and ovarian cancer panel  

1. Download the vcf file running from ahcg\_pipeline. [link](https://raw.githubusercontent.com/mjaeyi/ahcg_pipeline/master/variants.vcf)  

2. Download the golden vcf file running from GIAB.  
'$ wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/analysis/10XGenomics_calls_08142015/NA12878_phased_variants.vcf.gz'

3. Make the BED file for the genes from color genomics and otogenetics, breast and ovarian cancer panel.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Note: 1. [Gene list from otogenetics](http://www.otogenetics.com/forms/Breast\_Cancer\_gene\_list.pdf)    
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Note: 2. [Gene list from color genomics](https://s3.amazonaws.com/color-static-prod/pdfs/validationWhitePaper.pdf)    
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Note: 3. Manually make the gene list and look up their accession numbers in the reference genome file.    
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Note: [Ref Accession Number](http://www.genecards.org/cgi-bin/carddisp.pl?gene=PTEN)  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Note: The Ref Accession Number I am using is the first return after the "REFSEQ mRNAs" section.    
3.1. Make the gene\_name and ref\_accession\_number list    
&nbsp;&nbsp;&nbsp;&nbsp;Note: The first column is the gene symbols, the second column is the NCBI reference.     
&nbsp;&nbsp;&nbsp;&nbsp;Note: Name the output file as "BRC\_OC\_gene\_list.txt"    
3.2. Download the reference genome file    
`$ wget http://vannberg.biology.gatech.edu/data/ahcg2016/reference_genome/hg19_refGene.txt`     
3.3. Run the python pipeline to get the genome coordinates for the gene list  
&nbsp;&nbsp;&nbsp;&nbsp;Note: The python pipeline is in the folder GetCoordinatesForGenes.    
`$ python GetCoordinatesForGenes.py`  
&nbsp;&nbsp;&nbsp;&nbsp;Note: The example output is in the folder GetCoordinatesForGenes- "BRC\_OC\_gene\_list\_BED.txt".    
3.4. Add 20bp flank region to both the chrom\_star and chrom\_end in the "BRC\_OC\_gene\_list\_BED.txt" file.  
&nbsp;&nbsp;&nbsp;&nbsp;Note: The python pipeline is in the folder GetCoordinatesForGenes.    
`$ python Add 20bp.py`  
&nbsp;&nbsp;&nbsp;&nbsp;Note: The example output is in the folder GetCoordinatesForGenes- "BRC\_OC\_gene\_list\_BED\_add\_20.txt".  

4. Use bedtools to extract variants from ahcg\_pipeline vcf files and GIAB golden vcf, falling within the regions listed in the bed file for breast and ovarian cancer panel.   
`$ bedtools intersect -wa -header -a <vcf file> -b <bed file> > ahcg\_vcf\_BRCA\_OC.txt` 
`$ bedtools intersect -wa -header -a <vcf file> -b <bed file> > NA12878\_vcf\_BRCA\_OC.txt`  

5. Use bedtools to find intersect calls between vcf(from pipeline) and golden standard vcf.  
`$ bedtools intersect -a <vcf_golden_standard> -b <vcf_pipeline> > NA12878_ahcg_overlap_vcf_BRCA_OC.txt`   
  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Note: The output file could be found from vcf_files/NA12878_ahcg_overlap_vcf_BRCA_OC.txt   