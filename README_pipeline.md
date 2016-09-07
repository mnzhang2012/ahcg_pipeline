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
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;1.2.1. Reference genome
  Reference genomes can be downloaded from [Illumina iGenomes](http://support.illumina.com/sequencing/sequencing_software/igenome.html)  
 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;This pipeline is using [UCSC-hg19 reference genome]    
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`$ wget www.prism.gatech.edu/~sravishankar9/resources.tar.gz`

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;1.2.2. Test data
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
    The prebuilt index for reference genome hg19 could be download from [link](ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/hg19.zip)  

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
    &nbsp;&nbsp; 3. For the error of running out of memory for java application, you could change the give more memory for the virtual box and also make sure java's maximum memory usage has been increased according to how much you assigned for your virtual box using the arguments -Xmx<memory>.


###3. Help

To access help use the following command:
python3 ahcg_pipeline.py -h


[GATK - version 3.4](https://software.broadinstitute.org/gatk/download/)  
 [link](http://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/hg19.zip)  