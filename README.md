# Finding Dilated Cardiomyopathy(DCM) related genetic variants

####Project Goal: Finding Dilated Cardiomyopathy(DCM) related genetic variants      
Input: Two control and four DCM patients samples' bam file      
Output: Genetic variants with reads depth and pathogenic information.      
Procedues:     
1) Perfrom variant call analysis      
2) Read depth calculation on DCM-related gene panel       
3) Pathogenic annotation on DCM-related gene panel       
   
###1. Download the bam file for the 6 individuals.     
`$ sh get_DCM_sample.sh`    

###2. Perform variant call analysis using GATK-HaplotypeCaller.
` $ java -Xmx4g -jar lib/GenomeAnalysisTK.jar -T HaplotypeCaller -R input_files/ref_genome/hg19.fa -I input_files/dcm_bam/Control1_RG_MD_IR_BQ.bam --genotyping_mode DISCOVERY -stand_emit_conf 10 -stand_call_conf 30 -o input_files/dcm_vcf_haplotypeCaller/c1_raw_variants.vcf`    
###3. Recalibrated raw vcf using GATK-recalibrator.  
`$ java -Xmx4g -jar lib/GenomeAnalysisTK.jar -T VariantRecalibrator -R input_files/ref_genome/hg19.fa -input input_files/dcm_vcf_haplotypeCaller/c1_raw_variants.vcf -resource:hapmap,known=false,training=true,truth=true,prior=15.0 input_files/gatk_variantRecali/hapmap_3.3.hg19.sites.vcf -resource:omni,known=false,training=true,truth=false,prior=12.0 input_files/gatk_variantRecali/1000G_omni2.5.hg19.sites.vcf -resource:1000G,known=false,training=true,truth=false,prior=10.0 input_files/gatk_variantRecali/1000G_phase1.snps.high_confidence.hg19.sites.vcf -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 input_files/dbsnp/dbsnp_138.hg19.vcf -an DP -an QD -an FS -an SOR -an MQRankSum -an ReadPosRankSum -mode SNP -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 -recalFile input_files/dcm_vcf_haplotypeCaller/c1_recalibrate_SNP.recal -tranchesFile input_files/dcm_vcf_haplotypeCaller/c1_recalibrate_SNP.tranches`  

`$ java -Xmx4g -jar lib/GenomeAnalysisTK.jar -T ApplyRecalibration -R input_files/ref_genome/hg19.fa -input input_files/dcm_vcf_haplotypeCaller/c1_raw_variants.vcf -mode SNP --ts_filter_level 99.0 -recalFile input_files/dcm_vcf_haplotypeCaller/c1_recalibrate_SNP.recal -tranchesFile input_files/dcm_vcf_haplotypeCaller/c1_recalibrate_SNP.tranches -o c1_recalibrated_snps_raw_indels.vcf`  

&nbsp;&nbsp;Notes: Step 2 and Step 3 could be done in one step using script dcm_variant_call_pipeline.sh.     
`$ sh dcm_variant_call_pipeline.sh`    
###4. Create the bed file of the chromosome coordinates for DCM genes     
&nbsp;&nbsp;4.1 Create DCM related gene list    
&nbsp;&nbsp;&nbsp;&nbsp;Hint: Look up the DCM related publications and decide the genes that play important roles in DCM. The file has been saved as dcm\_gene.txt.   
&nbsp;&nbsp;4.2 Download the reference genome file    
`$ wget http://vannberg.biology.gatech.edu/data/ahcg2016/reference_genome/hg19_refGene.txt`     
&nbsp;&nbsp;4.3 Run the python script to get the genome coordinates for the NCBI\_Accession list      
`$ python get_coordinates_from_gene.py`    
&nbsp;&nbsp;4.4 Add 20bp flank region to both the chrom\_star and chrom\_end of each gene's chromosome coordinates.       
`$ python 20bp_flank_both_ends.py`     
&nbsp;&nbsp;&nbsp;&nbsp;Hint: The output file saved as dcm\_gene\_list.bed    

###5. Extract the variants of DCM genes for each patient samples    
`$ bedtools intersect -wa -header -a c1_recalibrated_snps_raw_indels.vcf -b dcm_gene_list.bed > control1_dcm_final.vcf`    

###6. Calculat the reads deapth information for DCM genes    
&nbsp;&nbsp;6.1 Extract brca1 alignments using samtools.    
`$ samtools view -L dcm_gene_list.bed Control1_RG_MD_IR_BQ.bam -b > c1_dcm.bam`     
&nbsp;&nbsp;6.2 Computes and summarize the depth for brca1    
`$ bedtools genomecov -ibam c1_dcm.bam -bga > c1_dcm_bga.bed`  
&nbsp;&nbsp;6.3 Intersection between two bed files.    
`$ bedtools intersect -split -a c1_dcm_bga.bed -b dcm_gene_list.bed -bed > c1_bcm.final.bed`  
&nbsp;&nbsp;6.4. Use the read\_depth\_calculation.py to get the read depth result.    
`$ python read_depth_calculation.py c1_bcm.final.bed > c1_bcm_depth.txt`   

&nbsp;&nbsp;&nbsp;&nbsp;Notes: Step 6 could be done for 6 samples using script dcm\_read\_coverage.sh.    
