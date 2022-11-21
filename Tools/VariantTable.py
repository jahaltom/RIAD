#Using output from Variant.py 
#wc -l *H* > genotypes

import pandas as pd

#Read in genotypes
gt=pd.read_csv("genotypes",sep='\s+',header=None)

vars=["rs1990760","rs60910145","rs73885319"]
for v in vars:
    #Loop through pops
    pop=["AFR","AMR","EUR","EAS","SAS"]
    table=[]
    for p in pop:
        #Variant of interest
        variant=gt[gt[1].str.contains(v)]
        #pop
        variant=variant[variant[1].str.contains(p)]
        #Get total(n) and calc genotype freq
        Total=variant[0].sum()
        variant[p] = variant[0] /Total   
        variant=variant.reset_index()     
        #(HomAlt*2 + Het)/Total*2. This will  be the alt allele freq.
        alt_alFreq=(variant[1:2][0][1]*2 + variant[0:1][0][0])/(Total*2)
        ref_alFreq=1-alt_alFreq
        ##Build table 
        #Genotype freq
        variant=variant.transpose()[p:p]
        variant.columns = ["Heterozygous","HomozygousAlt","HomozygousRef"]
        #Other columns
        variant["n"]=Total
        variant["Alt"]=alt_alFreq
        variant["Ref"]=ref_alFreq     
        table.append(variant)
    table=pd.concat(table)
    table['Super Population']=table.index
    table=table[["Super Population","n","Heterozygous","HomozygousAlt","HomozygousRef","Alt","Ref"]]
    table.to_csv(v+"_table.txt",sep='\t',index=False,mode='w')
