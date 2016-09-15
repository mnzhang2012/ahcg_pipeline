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