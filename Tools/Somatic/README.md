**TCGABiolinks.R:** Scans TCGA for studies with RNA-Seq Normal and Tumor with matching WES data. 
**GenerateMetadata.py:** Using output from TCGABiolinks.R, generates design.txt to be used with SomaticVarCall.py.
**RNA-SeqvsWES.table.py:** Downlaods COSMIC vcf. Makes a table compairing calls from SomaticVarCall.py and WES. 
**Boxplot.py:** Takes table from RNA-SeqvsWES.table.py.
