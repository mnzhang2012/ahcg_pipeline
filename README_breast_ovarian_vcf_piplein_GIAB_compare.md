# Compare the ahcg\_pipeline vcf file with GIAB vcf file for genes from color genomics and otogenetics, breast and ovarian cancer panel  

1. Download the vcf file running from ahcg\_pipeline. [link](https://raw.githubusercontent.com/mjaeyi/ahcg_pipeline/master/variants.vcf)  

2. Download the golden vcf file running from GIAB.  
'$ wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/analysis/10XGenomics_calls_08142015/NA12878_phased_variants.vcf.gz'

3. Make the BED file for the genes from color genomics and otogenetics, breast and ovarian cancer panel.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Note: 1. [Gene list from otogenetics](http://www.otogenetics.com/forms/Breast\_Cancer\_gene\_list.pdf)    
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Note: 2. [Gene list from color genomics](https://s3.amazonaws.com/color-static-prod/pdfs/validationWhitePaper.pdf)    
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Note: 3. Manually make the gene list and look up their accession numbers in the reference genome file.    
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Note: [Ref Accession Number](http://www.genecards.org/cgi-bin/carddisp.pl?gene=PTEN)  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Note: The Ref Accession Number I am using is the first return after the "REFSEQ mRNAs" section.    
3.1. Make the gene\_name and ref\_accession\_number list    
&nbsp;&nbsp;&nbsp;&nbsp;Note: The first column is the gene symbols, the second column is the NCBI reference.     
&nbsp;&nbsp;&nbsp;&nbsp;Note: Name the output file as "BRC\_OC\_gene\_list.txt"    
3.2. Download the reference genome file    
`$ wget http://vannberg.biology.gatech.edu/data/ahcg2016/reference_genome/hg19_refGene.txt`     
3.3. Run the python pipeline to get the genome coordinates for the gene list  
&nbsp;&nbsp;&nbsp;&nbsp;Note: The python pipeline is in the folder GetCoordinatesForGenes.    
`$ python GetCoordinatesForGenes.py`  
&nbsp;&nbsp;&nbsp;&nbsp;Note: The example output is in the folder GetCoordinatesForGenes- "BRC\_OC\_gene\_list\_BED.txt".    
3.4. Add 20bp flank region to both the chrom\_star and chrom\_end in the "BRC\_OC\_gene\_list\_BED.txt" file.  
&nbsp;&nbsp;&nbsp;&nbsp;Note: The python pipeline is in the folder GetCoordinatesForGenes.    
`$ python Add 20bp.py`  
&nbsp;&nbsp;&nbsp;&nbsp;Note: The example output is in the folder GetCoordinatesForGenes- "BRC\_OC\_gene\_list\_BED\_add\_20.txt".  

4. Use bedtools to extract variants from ahcg\_pipeline vcf files and GIAB golden vcf, falling within the regions listed in the bed file for breast and ovarian cancer panel.   
`$ bedtools intersect -wa -header -a <vcf file> -b <bed file> > ahcg\_vcf\_BRCA\_OC.txt` 
`$ bedtools intersect -wa -header -a <vcf file> -b <bed file> > NA12878\_vcf\_BRCA\_OC.txt`  

5. Use bedtools to find intersect calls between vcf(from pipeline) and golden standard vcf.  
`$ bedtools intersect -a <vcf_golden_standard> -b <vcf_pipeline> > NA12878_ahcg_overlap_vcf_BRCA_OC.txt`   
  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Note: The output file could be found from vcf_files/NA12878_ahcg_overlap_vcf_BRCA_OC.txt   