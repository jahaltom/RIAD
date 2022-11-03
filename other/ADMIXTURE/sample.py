import pandas as pd
from pandas import DataFrame


        
df=[]       
qfile=pd.read_csv("SAMPLE.prop",sep=' ',header=None)
qfile.columns =["EUR","EAS","AMR","SAS","AFR" ]############################Comfirm
prop=qfile.iloc[0].tolist()
max_prop_indx=prop.index(max(prop))
         
df.append("SAMPLE")
df.append(qfile.columns[max_prop_indx])    
df=pd.DataFrame(df).T
df.columns =["run_accession","Inferred Superpopulation"]
df=pd.concat([df, qfile], axis=1)
  ##############################################################################################################################Change header order match RIA       
df.to_csv('ADMIXTURE_Results.tsv',sep='\t',mode='w',index=False)     
df.to_csv('/work/LAS/xgu-lab/RIA/ADMIXTURE_Results.tsv',sep='\t',mode='a',index=False,header=None) 
