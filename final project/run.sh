#!/bin/bash
#!/usr/bin/Rscript
echo "Please use the following examples and then type in the name you are interested. The pipeline will generate the variant analysis" 

echo "The names of individual: Control1, Control2, Patient1, Patient2, Patient3, Patien4"
echo -n "Enter the name of individual > "
read text
echo "You entered: $text"
echo "The output variant file/analysis is the ${text}.final.vcf The reads depth coverage plot will be in the ${text} folder"

mkdir ${text}

echo "Genomic variant exploring pipeline starts!"
echo "Download bam file"

if [ ${text} = "Control1" ]
then
	wget http://vannberg.biology.gatech.edu/data/DCM/MenCo001DNA/Control1_RG_MD_IR_BQ.bam
	wget http://vannberg.biology.gatech.edu/data/DCM/MenCo001DNA/Control1_RG_MD_IR_BQ.bai
elif [ ${text} = "Control2" ]
then
	wget http://vannberg.biology.gatech.edu/data/DCM/MenCo002DNA/Control2_RG_MD_IR_BQ.bam
	wget http://vannberg.biology.gatech.edu/data/DCM/MenCo002DNA/Control2_RG_MD_IR_BQ.bai
elif [ ${text} = "Patient1" ]
then
	wget http://vannberg.biology.gatech.edu/data/DCM/MenPa001DNA/Patient1_RG_MD_IR_BQ.bam
	wget http://vannberg.biology.gatech.edu/data/DCM/MenPa001DNA/Patient1_RG_MD_IR_BQ.bai
elif [ ${text} = "Patient2" ]
then
	wget http://vannberg.biology.gatech.edu/data/DCM/MenPa002DNA/Patient2_RG_MD_IR_BQ.bam
	wget http://vannberg.biology.gatech.edu/data/DCM/MenPa002DNA/Patient2_RG_MD_IR_BQ.bai
elif [ ${text} = "Patient3" ]
then
	wget http://vannberg.biology.gatech.edu/data/DCM/MenPa003DNA/Patient3_RG_MD_IR_BQ.bam
	wget http://vannberg.biology.gatech.edu/data/DCM/MenPa003DNA/Patient3_RG_MD_IR_BQ.bai
elif [ ${text} = "Patient4" ]
then
	wget http://vannberg.biology.gatech.edu/data/DCM/MenPa004DNA/Patient4_RG_MD_IR_BQ.bam
	wget http://vannberg.biology.gatech.edu/data/DCM/MenPa004DNA/Patient4_RG_MD_IR_BQ.bai
else
	echo "Please try again"
	exit
fi

echo "Perform variant call analysis using GATK-HaplotyeCaller"
java -Xmx4g -jar lib/GenomeAnalysisTK.jar -T HaplotypeCaller -R input_files/ref_genome/hg19.fa -I ${text}_RG_MD_IR_BQ.bam --genotyping_mode DISCOVERY -stand_emit_conf 10 -stand_call_conf 30 -o ${text}/${text}_raw_variants.vcf
echo "Perform variants recalibration analysis using GATK-VariantRecalibrator"
java -Xmx4g -jar lib/GenomeAnalysisTK.jar -T VariantRecalibrator -R input_files/ref_genome/hg19.fa -input ${text}/${text}_raw_variants.vcf -resource:hapmap,known=false,training=true,truth=true,prior=15.0 input_files/gatk_variantRecali/hapmap_3.3.hg19.sites.vcf -resource:omni,known=false,training=true,truth=false,prior=12.0 input_files/gatk_variantRecali/1000G_omni2.5.hg19.sites.vcf -resource:1000G,known=false,training=true,truth=false,prior=10.0 input_files/gatk_variantRecali/1000G_phase1.snps.high_confidence.hg19.sites.vcf -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 input_files/dbsnp/dbsnp_138.hg19.vcf -an DP -an QD -an FS -an SOR -an MQRankSum -an ReadPosRankSum -mode SNP -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 -recalFile ${text}/${text}_recalibrate_SNP.recal -tranchesFile ${text}/${text}_recalibrate_SNP.tranches
java -Xmx4g -jar lib/GenomeAnalysisTK.jar -T ApplyRecalibration -R input_files/ref_genome/hg19.fa -input ${text}/${text}_raw_variants.vcf -mode SNP --ts_filter_level 99.0 -recalFile ${text}/${text}_recalibrate_SNP.recal -tranchesFile ${text}/${text}_recalibrate_SNP.tranches -o ${text}/${text}_recalibrated_snps_raw_indels.vcf

mkdir temp
echo "Create dcm-related 6 genes' genome position bed file"
./scripts/dcm_gene_bed.pl
echo "Computes and summarize the read depth coverage for dcm related genes "
samtools view -L dcm_exon_list.bed ${text}_RG_MD_IR_BQ.bam -b > ${text}/${text}.dcm.bam
echo "Reporting genome coverage"
bedtools genomecov -ibam ${text}/${text}.dcm.bam -bga > ${text}/${text}.dcm.bed
bedtools intersect -split -a ${text}/${text}.dcm.bed -b dcm_exon_list.bed -bed > ${text}/${text}.final.bed
python scripts/depth_for_each_site.py ${text}/${text}.final.bed ${text}/${text}.final.full.txt
echo "Calculateepth for each gene"
./scripts/label.pl ${text}/${text}.final.full.txt ${text}/${text}.final.full.labeled.txt

echo "Create reads coverage plot"
#sudo apt-get install r-base-core
#R
#install.packages('ggplot2' repos='http://cran.rstudio.com', type='source')
#q()

Rscript ./scripts/plots.R ${text}/${text}.final.full.labeled.txt

echo "Extract the dcm-related 6 genes' genomic variants"
bedtools intersect -wa -header -a ${text}/${text}_recalibrated_snps_raw_indels.vcf -b dcm_exon_list.bed > ${text}/${text}_dcm_final.vcf
echo "Annotate dcm_exon_list.bed with read depth and the exon number of each position and this will take a while..."
python scripts/position_gene_exon_depth.py dcm_gene_list.txt input_files/ref_genome/hg19_refGene.txt ${text}/${text}.final.full.txt > ${text}/${text}_dcm_position_exon_number_depth.txt
echo "Annotate individual's vcf file with the exon and depth information."
python scripts/add_infor_to_vcf.py ${text}/${text}_dcm_final.vcf ${text}/${text}_dcm_position_exon_number_depth.txt ${text}.final.vcf > ${text}/${text}.dcm.annotated.vcf
echo "${text}.final.vcf only included the variants(GATK-recalibrator passed) with depth greater than 20 for ${text}"
rm -r temp
rm cutoff*
rm gene_list_BED*
rm dcm_exon_list.bed
mv *.jpg ./${text}