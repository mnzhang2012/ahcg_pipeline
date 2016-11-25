# Variants calling pipeline for genomic data analysis

###1. Software Requirements  
1.1. Required programs   
    
  1. [Python3 - version 3.4.1](https://www.python.org/download/releases/3.4.1/)
  2. [Trimmomatic - version 0.36](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip)
  3. [Bowtie2 - version 2.2.9](https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.9/)
  4. [Picard tools - version 2.6.0](https://github.com/broadinstitute/picard/releases/download/2.6.0/picard.jar)
  5. [GATK - version 3.4](https://software.broadinstitute.org/gatk/download/)  
    
1.2. Requred Datasets   
&nbsp;&nbsp;1.2.1. Reference genome    
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Reference genomes(UCSC-hg19 reference genome) can be downloaded from [Illumina iGenomes](http://support.illumina.com/sequencing/sequencing_software/igenome.html)        
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;The genomic reference and dbSNP file needed for running pipeline could be download from [link](http://www.prism.gatech.edu/~sravishankar9/resources.tar.gz)    
    
&nbsp;&nbsp;1.2.2. Test dataset    
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Use the following commands to download and prepare test dataset(NIST sample NA12878)    

 >`$ wget ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R1_001.fastq.gz`  
  `$ wget ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/NIST7035_TAAGGCGA_L001_R2_001.fastq.gz`  
  `$ gunzip NIST7035_TAAGGCGA_L001_R1_001.fastq.gz`  
 ` $ gunzip NIST7035_TAAGGCGA_L001_R2_001.fastq.gz`  
  `$ head -1000 NIST7035_TAAGGCGA_L001_R1_001.fastq > test_r1.fastq`  
 ` $ head -1000 NIST7035_TAAGGCGA_L001_R2_001.fastq > test_r2.fastq`    

###2. Instructions for running pipeline    
2.1. Build bowite2 index for reference genome    
&nbsp;&nbsp;&nbsp;&nbsp;`$ bowtie2-build <reference_in> <bt2_base>`    
&nbsp;&nbsp;&nbsp;&nbsp;The prebuilt index for reference genome hg19 could be download from [link](ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/hg19.zip)    

2.2. Build fasta index for reference genome     
&nbsp;&nbsp;&nbsp;&nbsp;`$ samtools faidx ref.fasta`    

2.3. Build genome dict file for reference genome    
&nbsp;&nbsp;&nbsp;&nbsp;`$ java -jar picard.jar CreateSequenceDictonary R=ref.fasta O=ref.dict`
    
2.4. Reorder the VCF file   
&nbsp;&nbsp;&nbsp;&nbsp;`$ java -jar picard.jar SortVcf   I=original.vcf O=sorted.vcf SEQUENCE_DICTIONARY=ref.dict`    

2.5. Run ahcg_pipeline.py    
> `$ python3 ahcg_pipeline.py -t lib/Trimmomatic-0.36/trimmomatic-0.36.jar -b lib/bowtie2-2.2.9/bowtie2 -p lib/picard.jar -g lib/GenomeAnalysisTK.jar -i sequence/test_r1.fastq sequence/test_r2.fastq -w hg19 -d dbsnp/dbsnp_138.hg19.vcf -r genome/hg19.fa -a lib/Trimmomatic-0.36/adapters/NexteraPE-PE.fa -o output`     
&nbsp;&nbsp;&nbsp;&nbsp;Options:   
&nbsp;&nbsp;1. -i input files option  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; -i pair-end-read1 pair-end-read2  
&nbsp;&nbsp;2. -w bowtie2 index files option  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; -w hg19  

###3. Help  
To access help use the following command:  
`$ python3 ahcg_pipeline.py -h`    
    
    
    
# Change the remote url for git repositories
1. Change the current working directory to local git repository folder  
2. List exiting remotes for the repository    
   `$ git remote -v`
3. Use the git remote set-url command to change the remote repository's URL  
   `$ git remote set-url origin https://github.com/mnzhang2012/ahcg_pipeline.git`
4. Verify the remote URL has changed  
   `$ git remote -v`
    
    
# Update an existing file to a GitHub repository
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
    


# Extracting reads mapping to BRCA1 from NA12878 HiSeq Exome dataset  

1. Download the reference genome file  
`$ wget http://vannberg.biology.gatech.edu/data/ahcg2016/reference_genome/hg19_refGene.txt`  
  
2. Extract the genome coordinates for BRCA1  
`$ grep -w "NM_007294" hg19_refGene.txt > NM_007294.txt`  
&nbsp;&nbsp;Note: NM_007294 is the popular transcript of BRCA1  

3. Extract chrom, chromStart, chromEnd information from NM_007294 and convert it to bed file format.  
`$ awk '{print $3, $5-1000, $6+1000}' NM_007294.txt > NM_007294_gene_coordinates.bed` 

4. Download the NA12878 HiSeq Exome dataset BAM files from the [link](http://ftp://ftp-trace.ncbi.nih.gov/giab/ftp/data/NA12878/Garvan_NA12878_HG001_HiSeq_Exome/)    

5. Using samtools to subset the bam file to regions corresponding to BRCA1.  
> `$ samtools view -L <bed file> -b -o < outout bam file > < input bam file >`  
&nbsp;&nbsp;&nbsp;&nbsp;Note: -b just specifies that the output needs to be a bam file.   
&nbsp;&nbsp;&nbsp;&nbsp;Example code:
`$ samtools view -L NM_007294_gene_coordinates.bed -b -o BRCA1_1.bam project.NIST_NIST7035_H7AP8ADXX_TAAGGCGA_1_NA12878.bwa.markDuplicates.bam`  
`$ samtools view -L NM_007294_gene_coordinates.bed -b -o BRCA1_2.bam project.NIST_NIST7035_H7AP8ADXX_TAAGGCGA_2_NA12878.bwa.markDuplicates.bam`   
`$ samtools view -L NM_007294_gene_coordinates.bed -b -o BRCA1_3.bam project.NIST_NIST7086_H7AP8ADXX_TAAGGCGA_1_NA12878.bwa.markDuplicates.bam`   
`$ samtools view -L NM_007294_gene_coordinates.bed -b -o BRCA1_4.bam project.NIST_NIST7086_H7AP8ADXX_TAAGGCGA_2_NA12878.bwa.markDuplicates.bam`     

6. Merge bam files  
`$ samtools merge finalBamFile.bam *.bam`
  
7. Using bedtools to convert the bam file to a fastq file    
> `$ samtools sort -n aln.bam aln.qsort.bam`    
`$ bedtools bamtofastq -i <bam file> -fq < fastq r1> -fq2 < fastq r2>`   
Note:From the BRCA1 bam file we now extract the reads aligning to the region using bedtools      
Note:If your BAM alignments are from paired-end sequence data, one can use the -fq2 option to create two distinct FASTQ output files, one for end 1 and one for end 2. When using this option, it is required that the BAM file is sorted/grouped by the read name. This keeps the resulting records in the two output FASTQ files in the same order.   
8. Picard genome dictionary  
`$ java -jar picard.jar CreateSequenceDictionary R=ref.fa O=ref.dict`
&nbsp;&nbsp;Note: The reference file could be download from [link](http://vannberg.biology.gatech.edu/data/chr17.fa)    

9. Samtools fasta index   
`$ samtools faidx <path to ref.fa>`

10. Bowtie index  
`$ bowtie2-build <path to ref.fa> <output prefix>`

11. Run the ahcg\_pipeline.py to get BRCA1 variants    
`$ python3 ahcg_pipeline.py -t lib/Trimmomatic-0.36/trimmomatic-0.36.jar -b lib/bowtie2-2.2.9/bowtie2 -p lib/picard.jar -g lib/GenomeAnalysisTK.jar -i lab/BRCA1_bam/merge_BRCA1_end1.fq lab/BRCA1_bam/merge_BRCA1_end2.fq -w lab/chr17 -d dbsnp/dbsnp_138.hg19.vcf -r lab/chr17.fa -a lib/Trimmomatic-0.36/adapters/NexteraPE-PE.fa -o output_merge_BRC1`  

# Compare ahcg\_pipeline\_vcf with golden standard vcf file.   
1. Download the BED file, which contains the Nextera Rapid Capture Expanded Exome targeted regions' genomic coordinates.   
> `$ wget http://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/samplepreps_nextera/nexterarapidcapture/nexterarapidcapture_expandedexome_targetedregions.bed`  
  
2. Download the vcf file running from ahcg_pipeline.    
`$ wget http://vannberg.biology.gatech.edu/data/variants.vcf`    
  
3. Download the golden standard vcf file.  
> `$ wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/analysis/10XGenomics_calls_08142015/NA12878_phased_variants.vcf.gz`

4. Use bedtools to cut the Nextera Rapid Capture Expanded Exome targeted regions from golden standard vcf file and ahcg_pipeline vcf file.  
`$ bedtools intersect -wa -a <vcf file> -b <bed file>`  
  
5. Use bedtools to find intersect calls between ahcg_pipeline vcf and golden standard vcf.  
`$ bedtools intersect -header -a <vcf_golden_standard> -b <vcf_pipeline>`   
  
6. Use bedtools subtract command to find calls that are unique for gold_standard vcf.  
`$ bedtools substract -a <vcf_golden_standard> -b <vcf_intersect>`
  
7. Use bedtools subtract commandto find calls that are unique for ahcg_pipeline vcf.  
`$ bedtools substract -a <vcf_pipeline> -b <vcf_intersect>`
  
8. Use bedtools jaccard to find the number of TP, FP, FN variant calls.  
> `$ bedtools jaccard -a <vcf_intersect> -b <vcf_pipeline>`  
`$ bedtools jaccard -a <vcf_intersect> -b <vcf_golden_standard>`    

# Extract chromosome coordinates given NCBI_Accession number.    

1. Download the reference genome file  
`$ wget http://vannberg.biology.gatech.edu/data/ahcg2016/reference_genome/hg19_refGene.txt`   

2. Download the gene list of breast cancer biomarkers.  [link](http://www.otogenetics.com/forms/Breast_Cancer_gene_list.pdf)    
&nbsp;&nbsp;&nbsp;&nbsp;Note: Extract the first two columns and exported as a txt file.    
&nbsp;&nbsp;&nbsp;&nbsp;Note: The first column is the gene symbols, the second column is the NCBI reference.    
3. Run the python script to get the genome coordinates for the NCBI_Accession list    
`$ python get_coordinates_from_gene.py`  

# Compare ahcg\_pipeline\_vcf with GIAB\_vcf for BRC\_OC_genes

1. Download the vcf file running from ahcg\_pipeline. [link](https://raw.githubusercontent.com/mjaeyi/ahcg_pipeline/master/variants.vcf)  

2. Download the golden vcf file running from GIAB.  
`$ wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/latest/NA12878_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-Solid-10X_CHROM1-X_v3.3_highconf.vcf.gz`    
3. Make the BED file for the gene related with breast cancer(BRC) and ovarian cancer(OC)    
3.2. [Gene list from color genomics](https://www.otogenetics.com/forms/Breast_Cancer_gene_list.pdf)        
3.1. [Gene list from otogenetics](https://s3.amazonaws.com/color-static-prod/pdfs/validationWhitePaper.pdf)      
3.3. Create the gene list with gene symbols and associated transcript accession numbers.         
&nbsp;&nbsp;&nbsp;&nbsp;Note: The first column is the gene symbols, the second column is its associated transcript's NCBI reference number.     
3.4. Download the reference genome file    
`$ wget http://vannberg.biology.gatech.edu/data/ahcg2016/reference_genome/hg19_refGene.txt`     
3.5. Run the python script to get the genome coordinates for the gene list    
`$ python get_coordinates_from_gene.py`   
3.6. Add 20bp flank region to both the chrom\_star and chrom\_end of each gene's chromosome coordinates.       
`$ python 20bp_flank_both_ends.py`    
4. Use bedtools to extract variants from both ahcg\_pipeline vcf files and GIAB golden vcf, given the bed file of breast cancer and ovarian cancer related genes' chromosome coordinates.       
`$ bedtools intersect -wa -header -a <vcf file> -b <bed file> `       

5. Use bedtools to find intersect variants between vcf(from ahcg_pipeline) and golden standard vcf.  
`$ bedtools intersect -a <vcf_golden_standard> -b <vcf_ahcg_pipeline> `   


# How to apply GATK-variantRecalibrator on raw vcf file.  

1. Download the vcf file running from ahcg\_pipeline. [link](http://vannberglab.biology.gatech.edu/data/ahcg2016/vcf/NA12878_variants.vcf)      
  
2. Download the required known variant calls for VQSR to build model. [link](ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.8/hg19/)       
&nbsp;&nbsp;Note: Creating a tabix indexed vcf file    
> `$ gunzip file.gz`  
`$ bgzip -c file.vcf > file.vcf.gz`  
`$ tabix -p vcf file.vcf.gz`  
  
3. VariantRecalibrator - build a recalibration model to score variant quality for filtering purposes.    
> `$ java -Xmx4g -jar lib/GenomeAnalysisTK.jar  \  
    -T VariantRecalibrator \  
    -R ref_genome/hg19.fa  \  
    -input vcf_files_not_upload/NA12878_variants_ahcg.vcf \  
    -resource:hapmap,known=false,training=true,truth=true,prior=15.0 variantRecali/hapmap_3.3.hg19.sites.vcf \  -resource:omni,known=false,training=true,truth=false,prior=12.0 variantRecali/1000G_omni2.5.hg19.sites.vcf \  -resource:1000G,known=false,training=true,truth=false,prior=10.0 variantRecali/1000G_phase1.snps.high_confidence.hg19.sites.vcf \  -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 dbsnp/dbsnp_138.hg19.vcf 
    -an DP \  
    -an QD \  
    -an FS \  
    -an SOR \  
    -an MQRankSum \  
    -an ReadPosRankSum \  
    -mode SNP \  
    -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \  
    -recalFile recalibrate_SNP.recal \  
    -tranchesFile recalibrate_SNP.tranches`   
   
4. ApplyRecalibration - Apply a score cutoff to filter variants based on a recalibration table.    
> `$ java -jar lib/GenomeAnalysisTK.jar \ 
    -T ApplyRecalibration \ 
    -R ref_genome/hg19.fa \ 
    -input vcf_files_not_upload/NA12878_variants_ahcg.vcf \ 
    -mode SNP \ 
    --ts_filter_level 99.0 \ 
    -recalFile recalibrate_SNP.recal \ 
    -tranchesFile recalibrate_SNP.tranches \ 
    -o recalibrated_snps_raw_indels.vcf`
   

# How to calculate read depth based on alignment file(BRCA1)    
1. Extract BRCA1 gene chromosome coordinates from "brc\_oc\_gene\_list\_bed\_add\_20.txt"    
 `$ grep 'NM_007294' brc_oc_gene_list_bed_add_20.txt > brca1.bed`    
2. Extract brca1 alignments  
> `$ samtools view -L brca1.bed project.NIST_NIST7035_H7AP8ADXX_TAAGGCGA_1_NA12878.bwa.markDuplicates.bam -b > na12878.brca1.bam`  
       Note: -L: only output alignments overlapping in the input bed file   
       Note: -b: output alignments in the bam format   
  

3. Computes and summarize the depth for brca1  
> `$ bedtools genomecov -ibam na12878.brca1.bam -bga > na12878.brca1.bga.bed`   
       Note: -ibam: BAM file as input for coverage.   
       Note: -bga: Reporting genome coverage for all positions in BEDGRAPH format.    
       Note: This is BedGraph format:  chrom chromStart chromEnd dataValue  
4. Intersection between two bed files.    
> `$ bedtools intersect -split -a na12878.brca1.bga.bed -b brca1.bed -bed > brca1.final.bed`    
      Note: when using the -split option, only the exon overlaps are reported    
5. Use the read\_depth\_calculation.py to get the read depth result.    
`$ python read_depth_calculation.py`    

# How to annotate the vcf file with pathogenic information(BRCA1)         
1. Download the BRCA variant pathogenic annotation information file.    
> `$ wget http://vannberg.biology.gatech.edu/data/ahcg2016/BRCA/BRCA1_brca_exchange_variants.csv`  
`$ wget http://vannberg.biology.gatech.edu/data/ahcg2016/BRCA/BRCA2_brca_exchange_variants.csv`  

2. Use the pathgenic\_annotation\_add\_dp.py script to add pathogenic information to the variants.    
`$ python pathgenic_annotation_add_dp.py`  

# How to automate the pipeline to add read depth information and pathogenic information for the variants    
1. Run the ahcg\_pipeline with fastq file and create the vcf file.      
2. Run GATK genome recalibrator on the result vcf file.       
3. Use get\_coordinates\_from\_gene.py script to create the bed file for the chromosome coordinates of given genes of interests.    
4. Use bedtools to get the variants for the given bed file from the recalibrated vcf file.       
5. Use the pathgenic\_annotation\_add\_dp.py script to add pathogenic information to the variants.     
6. Use samtools to extract alignments given bed file and then use bedtools to compute and summarize the depth given bed file.    
7. Use read\_depth\_calculation.py script to calculate the reads depth based on bed file.   
5. Use the report\_depth\_pathogenic.py to integrate the reads depth information with the pathogenic information for the variants.    