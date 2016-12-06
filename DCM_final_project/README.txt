1. How to run the program
$ sh run.sh
Notes: The program will ask you to type in the name of individual you like to run.
Example of the name of individual: Control1, Control2, Patient1, Patient2, Patient3, Patient4

2. More information about the project and purpose of each step has shown below
1. General information
1) Project name: Investigate the dilated Cardiomyopathy(DCM) related genetic variants
2) Files
Input: Two control and four DCM patients samples' bam file      
Output: Genetic variants and reads depth coverage plot for each individual
3) Expectation   
1) Perfrom variant call analysis and variant recalibration     
2) Read depth calculation on DCM-related gene panel       
3) Reads depth coverage plot for DCM-related gene panel       

2. Procedues(Notes: All the codes below are using Patient1 as example)
1) Download the bam file for each of the 6 individuals.     
$ wget http://vannberg.biology.gatech.edu/data/DCM/MenPa001DNA/Patient1_RG_MD_IR_BQ.bam
$ wget http://vannberg.biology.gatech.edu/data/DCM/MenPa001DNA/Patient1_RG_MD_IR_BQ.bai      

2) Perform variant call analysis using GATK-HaplotypeCaller.
$ java -Xmx4g -jar lib/GenomeAnalysisTK.jar -T HaplotypeCaller -R input_files/ref_genome/hg19.fa -I Patient1_RG_MD_IR_BQ.bam --genotyping_mode DISCOVERY -stand_emit_conf 10 -stand_call_conf 30 -o Patient1/Patient1_raw_variants.vcf 

3) Recalibrated raw vcf using GATK-VariantRecalibrator
$ java -Xmx4g -jar lib/GenomeAnalysisTK.jar -T VariantRecalibrator -R input_files/ref_genome/hg19.fa -input Patient1/Patient1_raw_variants.vcf -resource:hapmap,known=false,training=true,truth=true,prior=15.0 input_files/gatk_variantRecali/hapmap_3.3.hg19.sites.vcf -resource:omni,known=false,training=true,truth=false,prior=12.0 input_files/gatk_variantRecali/1000G_omni2.5.hg19.sites.vcf -resource:1000G,known=false,training=true,truth=false,prior=10.0 input_files/gatk_variantRecali/1000G_phase1.snps.high_confidence.hg19.sites.vcf -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 input_files/dbsnp/dbsnp_138.hg19.vcf -an DP -an QD -an FS -an SOR -an MQRankSum -an ReadPosRankSum -mode SNP -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 -recalFile Patient1/Patient1_recalibrate_SNP.recal -tranchesFile Patient1/Patient1_recalibrate_SNP.tranches

4) Apply VQSR model using GATK-ApplyRecalibration
$ java -Xmx4g -jar lib/GenomeAnalysisTK.jar -T ApplyRecalibration -R input_files/ref_genome/hg19.fa -input Patient1/Patient1_raw_variants.vcf -mode SNP --ts_filter_level 99.0 -recalFile Patient1/Patient1_recalibrate_SNP.recal -tranchesFile Patient1/Patient1_recalibrate_SNP.tranches -o Patient1/Patient1_recalibrated_snps_raw_indels.vcf

5) Create the chromosomal coordinates bed file for DCM genes     
$ perl ./scripts/dcm_gene_bed.pl
The bed file for DCM gene is names "dcm_exon_list.bed"

6) Calculat the reads depth information for DCM genes and create the reads depth plot   
6.1) Extract alignments of dcm related genes using samtools
$ samtools view -L dcm_exon_list.bed Patient1_RG_MD_IR_BQ.bam -b > Patient1/Patient1.dcm.bam
6.2) Computes and summarize the read depth coverage for dcm related genes   
$ bedtools genomecov -ibam Patient1/Patient1.dcm.bam -bga > Patient1/Patient1.dcm.bed
6.3) Read coverage calculation for each gene 
$ bedtools intersect -split -a Patient1/Patient1.dcm.bed -b dcm_exon_list.bed -bed > Patient1/Patient1.final.bed
$ python scripts/depth_for_each_site.py Patient1/Patient1.final.bed Patient1/Patient1.final.full.txt
$ perl ./scripts/label.pl Patient1/Patient1.final.full.txt Patient1/Patient1.final.full.labeled.txt

7) Coverage plot generation
7.1) Install R and package ggplot2
sudo apt-get install r-base-core
R
install.packages('ggplot2' repos='http://cran.rstudio.com', type='source')
q()
7.2) Run plots.R to generate coverage plot for each gene
$ Rscript ./scripts/plots.R Patient1/Patient1.final.full.labeled.txt

8) Dcm related vcf analysis based GATK-variantRecalibrator and read_depth information
8.1) Extract the variants of DCM genes for each patient samples giving bed file    
$ bedtools intersect -wa -header -a Patient1/Patient1_recalibrated_snps_raw_indels.vcf -b dcm_exon_list.bed > Patient1/Patient1_dcm_final.vcf
8.2) Annotate dcm_exon_list.bed with read depth and the exon number of each position
$ python scripts/position_gene_exon_depth.py dcm_gene_list.txt input_files/ref_genome/hg19_refGene.txt Patient1/Patient1.final.full.txt > Patient1/Patient1_dcm_position_exon_number_depth.txt
8.3) Annotate individual's vcf file with the exon and depth information
$ python scripts/add_infor_to_vcf.py Patient1/Patient1_dcm_final.vcf Patient1/Patient1_dcm_position_exon_number_depth.txt Patient1.final.vcf > Patient1/Patient1.dcm.annotated.vcf
The final vcf for patient1 only included the variants(GATK-recalibrator passed) with depth greater than 20.