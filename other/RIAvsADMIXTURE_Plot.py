import pandas as pd
import matplotlib.pyplot as plt
import os


os.chdir(r"C:\Users\15154\Documents")

#Read in runtimes and file sizes.
RIA=pd.read_csv('NoNoiseRIA.tsv',sep='\t',header=None)  
RIA.columns =['SampleID', 'Runtime']
ADMIXTURE=pd.read_csv('NoNoiseADMIX.tsv',sep='\t',header=None)  
ADMIXTURE.columns =['SampleID', 'Runtime']
sizes=pd.read_csv('sizes.tsv',sep='\t',header=None)  
sizes.columns =['SampleID', 'SizeMB']

#Merge runtimes and file sizes
RIA=pd.merge(sizes,RIA,on=["SampleID"])
ADMIXTURE=pd.merge(sizes,ADMIXTURE,on=["SampleID"])



#Plot

plt.scatter("SizeMB","Runtime",data=RIA,s=5)
plt.scatter("SizeMB","Runtime",data=ADMIXTURE,s=5)
plt.legend(["RIA", "ADMIXTURE"], loc=0)
plt.ylabel("Runtime (Hours) ")
plt.xlabel("Compressed VCF File Size (MB) ")


plt.savefig("RIAvsADMIXTURE_RuntimeFileSize.png",format='png',dpi=150,bbox_inches='tight')
