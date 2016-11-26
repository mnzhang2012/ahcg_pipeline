#!/bin/bash

samtools view -L dcm_gene_list.bed Control1_RG_MD_IR_BQ.bam -b > c1_dcm.bam
bedtools genomecov -ibam c1_dcm.bam -bga > c1_dcm_bga.bed
bedtools intersect -split -a c1_dcm_bga.bed -b dcm_gene_list.bed -bed > c1_dcm.final.bed
python read_depth_calculation.py c1_dcm.final.bed > c1_dcm_depth.txt

samtools view -L dcm_gene_list.bed Control2_RG_MD_IR_BQ.bam -b > c2_dcm.bam
bedtools genomecov -ibam c2_dcm.bam -bga > c2_dcm_bga.bed
bedtools intersect -split -a c2_dcm_bga.bed -b dcm_gene_list.bed -bed > c2_dcm.final.bed
python read_depth_calculation.py c2_dcm.final.bed > c2_dcm_depth.txt

samtools view -L dcm_gene_list.bed Patient1_RG_MD_IR_BQ.bam -b > p1_dcm.bam
bedtools genomecov -ibam p1_dcm.bam -bga > p1_dcm_bga.bed
bedtools intersect -split -a p1_dcm_bga.bed -b dcm_gene_list.bed -bed > p1_dcm.final.bed
python read_depth_calculation.py p1_dcm.final.bed > p1_dcm_depth.txt

samtools view -L dcm_gene_list.bed Patient2_RG_MD_IR_BQ.bam -b > p2_dcm.bam
bedtools genomecov -ibam p2_dcm.bam -bga > p2_dcm_bga.bed
bedtools intersect -split -a p2_dcm_bga.bed -b dcm_gene_list.bed -bed > p2_dcm.final.bed
python read_depth_calculation.py p2_dcm.final.bed > p2_dcm_depth.txt

samtools view -L dcm_gene_list.bed Patient3_RG_MD_IR_BQ.bam -b > p3_dcm.bam
bedtools genomecov -ibam p3_dcm.bam -bga > p3_dcm_bga.bed
bedtools intersect -split -a p3_dcm_bga.bed -b dcm_gene_list.bed -bed > p3_dcm.final.bed
python read_depth_calculation.py p3_dcm.final.bed > p3_dcm_depth.txt

samtools view -L dcm_gene_list.bed Patient4_RG_MD_IR_BQ.bam -b > p4_dcm.bam
bedtools genomecov -ibam p4_dcm.bam -bga > p4_dcm_bga.bed
bedtools intersect -split -a p4_dcm_bga.bed -b dcm_gene_list.bed -bed > p4_dcm.final.bed
python read_depth_calculation.py p4_dcm.final.bed > p4_dcm_depth.txt