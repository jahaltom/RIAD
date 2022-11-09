**Variant.py:** Given a python dictionary of varaints as below, this script will search through genomic data and report the genotype of the sample at the variant site. 
Searches an dreports by self-reported ethnicity. Also needs metadata.txt containing self-reported ethnicities.

rs73885319 and rs60910145 https://clinvarminer.genetics.utah.edu/variants-by-gene/APOL1/significance/risk%20factor

rs1990760:  https://www.ncbi.nlm.nih.gov/snp/rs1990760


```
vars={
        #VariantID, Chr, Position, Ref, Alt.
        "rs1990760":	[2, 162267541 ,'C','T'],
        "rs73885319":	[22, 36265860,'A','G'],
        "rs60910145":	[22, 36265988,'T','G']
      
}
```

**POP:** a file like below.

```
SAS
AFR
AMR
EUR
EAS
```

To run use the folowing:

```
cat POP | while read i;do
 	sed "s/POP/$i/g" snakefile  > snakefile$i 
  snakemake -j 300 -k -s snakefile$i --cluster "sbatch -t 01:00:00 -c 10 -p RM-shared"
done
```
