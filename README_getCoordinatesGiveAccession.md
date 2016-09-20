# Extracting chromosome coordinates given NCBI_Accession number.    

1. Download the reference genome file  
`$ wget http://vannberg.biology.gatech.edu/data/ahcg2016/reference_genome/hg19_refGene.txt`   

2. Download the list of breast cancer biomarkers.  
`$ wget http://www.otogenetics.com/forms/Breast_Cancer_gene_list.pdf`   
&nbsp;&nbsp;&nbsp;&nbsp;Note: Extract the first two columns and exported as a txt file.    
&nbsp;&nbsp;&nbsp;&nbsp;Note: The first column is the gene symbols, the second column is the NCBI reference.    

3. Run the python pipeline  
&nbsp;&nbsp;&nbsp;&nbsp;Note: The python pipeline is in the folder GetCoordinatesForGenes.  
`$ python GetCoordinatesForGenes.py`  
&nbsp;&nbsp;&nbsp;&nbsp;Note: The example output is in the folder GetCoordinatesForGenes.  