#!/bin/bash

wget http://vannberg.biology.gatech.edu/data/DCM/MenCo001DNA/Control1_RG_MD_IR_BQ.bam
wget http://vannberg.biology.gatech.edu/data/DCM/MenCo001DNA/Control1_RG_MD_IR_BQ.bai

wget http://vannberg.biology.gatech.edu/data/DCM/MenCo002DNA/Control2_RG_MD_IR_BQ.bam
wget http://vannberg.biology.gatech.edu/data/DCM/MenCo002DNA/Control2_RG_MD_IR_BQ.bai

wget http://vannberg.biology.gatech.edu/data/DCM/MenPa001DNA/Patient1_RG_MD_IR_BQ.bam
wget http://vannberg.biology.gatech.edu/data/DCM/MenPa001DNA/Patient1_RG_MD_IR_BQ.bai

wget http://vannberg.biology.gatech.edu/data/DCM/MenPa002DNA/Patient2_RG_MD_IR_BQ.bam
wget http://vannberg.biology.gatech.edu/data/DCM/MenPa002DNA/Patient2_RG_MD_IR_BQ.bai

wget http://vannberg.biology.gatech.edu/data/DCM/MenPa003DNA/Patient3_RG_MD_IR_BQ.bam
wget http://vannberg.biology.gatech.edu/data/DCM/MenPa003DNA/Patient3_RG_MD_IR_BQ.bai

wget http://vannberg.biology.gatech.edu/data/DCM/MenPa004DNA/Patient4_RG_MD_IR_BQ.bam
wget http://vannberg.biology.gatech.edu/data/DCM/MenPa004DNA/Patient4_RG_MD_IR_BQ.bai


md5sum Control1_RG_MD_IR_BQ.bam Control1_RG_MD_IR_BQ.bai Control2_RG_MD_IR_BQ.bam Control2_RG_MD_IR_BQ.bai Patient1_RG_MD_IR_BQ.bam Patient1_RG_MD_IR_BQ.bai Patient2_RG_MD_IR_BQ.bam Patient2_RG_MD_IR_BQ.bai Patient3_RG_MD_IR_BQ.bam Patient3_RG_MD_IR_BQ.bai Patient4_RG_MD_IR_BQ.bam Patient4_RG_MD_IR_BQ.bai > md5sums.md5
md5sum -c md5sums.md5


