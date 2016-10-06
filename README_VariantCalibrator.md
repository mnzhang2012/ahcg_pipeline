# How to use GATK-variantRecalibrator to do vcf quality control     

1. Download the vcf file running from ahcg\_pipeline(reanalysis by using fastq file). [link](http://vannberglab.biology.gatech.edu/data/ahcg2016/vcf/NA12878_variants.vcf)      
  
2. Download the required known variatn calls for VQSR to build model. [link](ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.8/hg19/)       
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Note: Creating a tabix indexed vcf file    
`$ gunzip file.gz`  
`$ bgzip -c file.vcf > file.vcf.gz`  
`$ tabix -p vcf file.vcf.gz`  
  
3. Use GATK variant recalibrator on vcf file.    
`$ java -Xmx4g -jar lib/GenomeAnalysisTK.jar -T VariantRecalibrator -R ref_genome/hg19.fa -input vcf_files_not_upload/NA12878_variants_ahcg.vcf -resource:hapmap,known=false,training=true,truth=true,prior=15.0 variantRecali/hapmap_3.3.hg19.sites.vcf -resource:omni,known=false,training=true,truth=false,prior=12.0 variantRecali/1000G_omni2.5.hg19.sites.vcf -resource:1000G,known=false,training=true,truth=false,prior=10.0 variantRecali/1000G_phase1.snps.high_confidence.hg19.sites.vcf -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 dbsnp/dbsnp_138.hg19.vcf -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an InbreedingCoeff -mode SNP -recalFile output.recal -tranchesFile output.tranches -rscriptFile output.plots.R`   

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Note: "-an InbreedingCoeff" removed. Running the command with this tag will bring the error:  "Bad input: Values for InbreedingCoeff annotation not detected for ANY training variant in the input callset. VariantAnnotator may be usep://gatkforums.broadinstitute.org/discussion/49/using-variant-annotator"  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Note: The output files are in the folder "variantRecali-output"  