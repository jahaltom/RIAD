import pandas as pd
from pandas import DataFrame




#Read in metadata
metadata=pd.read_csv("metadata",sep="\t")
#Read in SNPs
snps=pd.read_csv("AllRunsSNPS.txt",sep="\t",header=None)   
snps.columns =['run_accession','Num of SNPs']
#Read in Unique mapped read Num
readNum=pd.read_csv("AllRunsNumUniqReads.txt",sep="\t",header=None)  
readNum.columns =['run_accession','Num of Unique Mapped Reads']
#Read in ancestry inference results 
results=pd.read_csv("SuperpopulationChrAll.PC20SVMResults",sep="\t")       

#Merge metadata with results  
resultsMeta=pd.merge(results,metadata,on=['run_accession'])  

#Merge with snps  
resultsMeta=pd.merge(resultsMeta,snps,on=['run_accession'])  
     
#Merge with Unique mapped read Num
resultsMeta=pd.merge(resultsMeta,readNum,on=['run_accession'])  
              

#Set whole field to No then change based off accuracy 
resultsMeta['Accurate'] ="No"
for idn in resultsMeta.index:
    if ((resultsMeta['Superpopulation'][idn] == resultsMeta['Eth1'][idn]) | (resultsMeta['Superpopulation'][idn] == resultsMeta['Eth2'][idn]))==True:
        resultsMeta['Accurate'][idn] ="Yes"

    
    
    

resultsMeta.to_csv('SNPs_ReadNum_metadata.tsv',sep='\t',mode='w',index=False) 

