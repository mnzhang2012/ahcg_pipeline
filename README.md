# Variant calling pipeline for genomic data analysis

###1. Software Requirements  
1.1. Here is a list of the programs that are needed for running this pipeline  
   &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Download the programs using wget command

  1. [Python3 - version 3.4.1](https://www.python.org/download/releases/3.4.1/)
  2. [Trimmomatic - version 0.36](http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip)
  3. [Bowtie2 - version 2.2.9](https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.9/)
  4. [Picard tools - version 2.6.0](https://github.com/broadinstitute/picard/releases/download/2.6.0/picard.jar)
  5. [GATK - version 3.4](https://software.broadinstitute.org/gatk/download/)  

1.2. Reference genome and Test data  
 &nbsp;&nbsp;1.2.1. Reference genome
  Reference genomes can be downloaded from [Illumina iGenomes](http://support.illumina.com/sequencing/sequencing_software/igenome.html)  
 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;This pipeline is using [UCSC-hg19 reference genome]    
 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;The genomic reference and dbSNP file that you need for the pipeline are avaiable at following link
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;`$ wget www.prism.gatech.edu/~sravishankar9/resources.tar.gz`

 &nbsp;&nbsp;1.2.2. Test data  
 &nbsp;&nbsp; Use the following protocol to download and prepare test dataset from NIST sample NA12878  

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

2.5. Run the ahcg_pipeline.py script  
   &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Exmaple command  
   > `$ python3 ahcg_pipeline.py -t lib/Trimmomatic-0.36/trimmomatic-0.36.jar -b lib/bowtie2-2.2.9/bowtie2 -p lib/picard.jar -g lib/GenomeAnalysisTK.jar -i sequence/test_r1.fastq sequence/test_r2.fastq -w hg19 -d dbsnp/dbsnp_138.hg19.vcf -r genome/hg19.fa -a lib/Trimmomatic-0.36/adapters/NexteraPE-PE.fa -o output`   
   &nbsp;&nbsp;Notes:   
     &nbsp;&nbsp;1. -i input files option  
       &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; -i pair-end-read1 pair-end-read2  
     &nbsp;&nbsp;2. -w bowtie2 index files option  
       &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; -w hg19  

###3. Help

To access help use the following command:
python3 ahcg_pipeline.py -h



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
Notes: origin is an alias on your system for a particular remote repository. It's not actually a property of that repository.    
Notes: A branch in Git is simply a lightweight movable pointer to one of these commits. The default branch name in Git is master. As you initially make commits, you're given a master branch that points to the last commit you made. Every time you commit, it moves forward automatically.
  
  
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
`$ samtools view -L <bed file> -b -o < outout bam file > < input bam file >`  
&nbsp;&nbsp;&nbsp;&nbsp;Note: -b just specifies that the output needs to be a bam file.   
> Example code:
`$ samtools view -L NM_007294_gene_coordinates.bed -b -o BRCA1_1.bam project.NIST_NIST7035_H7AP8ADXX_TAAGGCGA_1_NA12878.bwa.markDuplicates.bam`  
`$ samtools view -L NM_007294_gene_coordinates.bed -b -o BRCA1_2.bam project.NIST_NIST7035_H7AP8ADXX_TAAGGCGA_2_NA12878.bwa.markDuplicates.bam`   
`$ samtools view -L NM_007294_gene_coordinates.bed -b -o BRCA1_3.bam project.NIST_NIST7086_H7AP8ADXX_TAAGGCGA_1_NA12878.bwa.markDuplicates.bam`   
`$ samtools view -L NM_007294_gene_coordinates.bed -b -o BRCA1_4.bam project.NIST_NIST7086_H7AP8ADXX_TAAGGCGA_2_NA12878.bwa.markDuplicates.bam` 
6. Merge bam files  
`$ samtools merge finalBamFile.bam *.bam`
  
7. Using bedtools to convert the bam file to a fastq file    
`$ samtools sort -n aln.bam aln.qsort.bam`    
`$ bedtools bamtofastq -i <bam file> -fq < fastq r1> -fq2 < fastq r2>`  
&nbsp;&nbsp;Note:From the BRCA1 bam file we now extract the reads aligning to the region using bedtools      
&nbsp;&nbsp;Note:If your BAM alignments are from paired-end sequence data, one can use the -fq2 option to create two distinct FASTQ output files, one for end 1 and one for end 2. When using this option, it is required that the BAM file is sorted/grouped by the read name. This keeps the resulting records in the two output FASTQ files in the same order.   
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

# Compare with golden standard vcf file.   
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

# Extract chromosome coordinates given NCBI_Accession number.    

1. Download the reference genome file  
`$ wget http://vannberg.biology.gatech.edu/data/ahcg2016/reference_genome/hg19_refGene.txt`   

2. Download the list of breast cancer biomarkers.  
`$ wget http://www.otogenetics.com/forms/Breast_Cancer_gene_list.pdf`   
&nbsp;&nbsp;&nbsp;&nbsp;Note: Extract the first two columns and exported as a txt file.    
&nbsp;&nbsp;&nbsp;&nbsp;Note: The first column is the gene symbols, the second column is the NCBI reference.    
3. Run the python pipeline  
`$ python GetCoordinatesForGenes.py`  

# Compare with GIAB\_vcf\_file for BRC\_OC_genes

1. Download the vcf file running from ahcg\_pipeline. [link](https://raw.githubusercontent.com/mjaeyi/ahcg_pipeline/master/variants.vcf)  

2. Download the golden vcf file running from GIAB.  
`$ wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/latest/NA12878_GIAB_highconf_CG-IllFB-IllGATKHC-Ion-Solid-10X_CHROM1-X_v3.3_highconf.vcf.gz`

3. Make the BED file for the genes related with BRC and OC.  
3.1. [Gene list from otogenetics](http://www.otogenetics.com/forms/Breast\_Cancer\_gene\_list.pdf)    
3.2. [Gene list from color genomics](https://s3.amazonaws.com/color-static-prod/pdfs/validationWhitePaper.pdf)    
3.3. Manually make the gene list and look up their accession numbers in the reference genome file.    
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Note: [Ref Accession Number](http://www.genecards.org/cgi-bin/carddisp.pl?gene=PTEN)  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Note: The Ref Accession Number I am using is the first return after the "REFSEQ mRNAs" section.    
3.4. Make the gene\_name and ref\_accession\_number list    
&nbsp;&nbsp;&nbsp;&nbsp;Note: The first column is the gene symbols, the second column is the NCBI reference.     
&nbsp;&nbsp;&nbsp;&nbsp;Note: Name the output file as "BRC\_OC\_gene\_list.txt"    
3.5. Download the reference genome file    
`$ wget http://vannberg.biology.gatech.edu/data/ahcg2016/reference_genome/hg19_refGene.txt`     
3.6. Run the python script to get the genome coordinates for the gene list    
`$ python GetCoordinatesForGenes.py`   
3.7. Add 20bp flank region to both the chrom\_star and chrom\_end.  
`$ python Add 20bp.py`    
4. Use bedtools to extract variants from ahcg\_pipeline vcf files and GIAB golden vcf, given the regions listed in the bed file for breast and ovarian cancer panel.   
`$ bedtools intersect -wa -header -a <vcf file> -b <bed file> `   
`$ bedtools intersect -wa -header -a <vcf file> -b <bed file> `  

5. Use bedtools to find intersect calls between vcf(from pipeline) and golden standard vcf.  
`$ bedtools intersect -a <vcf_golden_standard> -b <vcf_pipeline> `   


# GATK-variantRecalibrator

1. Download the vcf file running from ahcg\_pipeline. [link](http://vannberglab.biology.gatech.edu/data/ahcg2016/vcf/NA12878_variants.vcf)      
  
2. Download the required known variant calls for VQSR to build model. [link](ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.8/hg19/)       
&nbsp;&nbsp;Note: Creating a tabix indexed vcf file    
`$ gunzip file.gz`  
`$ bgzip -c file.vcf > file.vcf.gz`  
`$ tabix -p vcf file.vcf.gz`  
  
3. Run the variant recalibrator.    
`$ java -Xmx4g -jar lib/GenomeAnalysisTK.jar  \  
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
   
4. Run GATK again to get vcf file.    
`$ java -jar lib/GenomeAnalysisTK.jar \ 
    -T ApplyRecalibration \ 
    -R ref_genome/hg19.fa \ 
    -input vcf_files_not_upload/NA12878_variants_ahcg.vcf \ 
    -mode SNP \ 
    --ts_filter_level 99.0 \ 
    -recalFile recalibrate_SNP.recal \ 
    -tranchesFile recalibrate_SNP.tranches \ 
    -o recalibrated_snps_raw_indels.vcf`
   

# How to calculate read depth based on alignment file.
 1. Extract BRCA1 gene chromosome coordinates from "BRC\_OC\_gene\_list\_BED.txt"  
 `$ grep 'NM_007298' brc_oc_gene_list_bed.txt > brca1.bed`
 2. Extract brca1 alignments  
 `$ samtools view -L brca1.bed project.NIST_NIST7035_H7AP8ADXX_TAAGGCGA_1_NA12878.bwa.markDuplicates.bam -b > na12878.brca1.bam`   
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Note: -L: only output alignments overlapping in the input bed file   
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Note: -b: output alignments in the bam format   
 3. Computes and summarize coverage for brca1   
 `$ bedtools genomecov -ibam na12878.brca1.bam -bga > na12878.brca1.bga.bed`   
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Note: -ibam BAM file as input for coverage.   
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Note: -bga  Reporting genome coverage for all positions in BEDGRAPH format.    
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Note: This is BedGraph format:  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;chrom chromStart chromEnd dataValue  
 4. Intersection between two bed files.  
`$ bedtools intersect -split -a na12878.brca1.bga.bed -b brca1.bed -bed > brca1.final.bed`    
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Note: when using the -split option, only the exon overlaps are reported    
 5. Use the ReadDepthCalculation.py to get the read depth result.    
`$ python ReadDepthCalculation.py`    

# How to annotate the vcf file with pathogenicity     
1. Download the BRCA variant pathogenic annotation file.    
`$ wget http://vannberg.biology.gatech.edu/data/ahcg2016/BRCA/BRCA1_brca_exchange_variants.csv`  
`$ wget http://vannberg.biology.gatech.edu/data/ahcg2016/BRCA/BRCA2_brca_exchange_variants.csv`  

2. Use the PathgenicAnnotation.py script to match the variants between vcf file and pathogenic annotation file.    
`$ python PathgenicAnnotation.py variantsRecali_output/recalibrated_snps_raw_indels.vcf variantsRecali_output/BRCA1_brca_exchange_variants.csv > variants_annotation.vcf`  

# How to automate the pipeline from given gene name to pathogenic variants   
1. Run the pipeline with fastq file and then got the vcf file.  
2. Run GATK genome recalibrator on vcf.  
3. Match vcf file with clinical variants information in order to find the deleterious variants.  
4. Calculate the reads coverage in order to prove the occurence of negative results is not due to low reads coverage. 
5. Integrate the reads coverage information with the pathogenic information for each variants.   
`$ python brca1_report_depth_pathogenic.py variants_annotation.vcf brca1_depth.txt`
