import pandas as pd
import glob
import matplotlib.pyplot as plt
import numpy as np




files=glob.glob("*Total*")

table = pd.DataFrame()
for i in files:
    df=pd.read_csv(i,sep='\t',header=None)
    df.columns =[i]
    table=pd.concat([table, df], axis=1)
        
        
files=glob.glob("*RNA_DNA*")

for i in files:
    df=pd.read_csv(i,sep='\t',header=None)
    df.columns =[i]
    table=pd.concat([table, df], axis=1)
        
     
        
table=table[[ 'TotalRNACalls','TotalDNACalls','TotalRNA_PASS_Calls' , 'TotalDNA_PASS_Calls', 'TotalRNA.COSMIC.Calls', 'TotalDNA.COSMIC.Calls','TotalRNA.COSMIC._PASS_Calls','TotalDNA.COSMIC._PASS_Calls','RNA_DNA_Calls','RNA_DNA_PASS_Calls']]



table=table[["TotalRNACalls","TotalDNACalls"]]
table = table.set_axis(["RNA-Seq Unfiltered","WES Unfiltered"], axis=1, inplace=False)

vals, names, xs = [],[],[]
for i, col in enumerate(table.columns):
    vals.append(table[col].values)
    names.append(col)
    xs.append(np.random.normal(i + 1, 0.04, table[col].values.shape[0]))  # adds jitter to the data points - can be adjusted
    
plt.boxplot(vals, labels=names)
palette = ['b', 'g', 'b', 'g','b', 'g', 'b', 'g','b', 'g']
for x, val, c in zip(xs, vals, palette):
    plt.scatter(x, val, alpha=0.4, color=c)
plt.xticks(rotation=45, ha='right')
plt.ylabel ('# Variant Calls')
plt.savefig("Total_RNA_WES_Calls_Unfiltered.png",format='png',dpi=150,bbox_inches='tight')



table=pd.read_csv("Somatic.txt",sep='\t')
table=table[['TotalRNA_PASS_Calls' , 'TotalDNA_PASS_Calls', 'TotalRNA.COSMIC.Calls', 'TotalDNA.COSMIC.Calls','TotalRNA.COSMIC._PASS_Calls','TotalDNA.COSMIC._PASS_Calls','RNA_DNA_Calls','RNA_DNA_PASS_Calls']]
table = table.set_axis(["RNA-Seq PASS","WES PASS","RNA-Seq COSMIC","WES COSMIC","RNA-Seq COSMIC PASS","WES COSMIC PASS","Common RNA-Seq/WES Unfiltered","Common RNA-Seq/WES PASS"], axis=1, inplace=False)
vals, names, xs = [],[],[]
for i, col in enumerate(table.columns):
    vals.append(table[col].values)
    names.append(col)
    xs.append(np.random.normal(i + 1, 0.04, table[col].values.shape[0]))  # adds jitter to the data points - can be adjusted
    
plt.boxplot(vals, labels=names)
palette = ['b', 'g', 'b', 'g','b', 'g', 'b', 'g','b', 'g']
for x, val, c in zip(xs, vals, palette):
    plt.scatter(x, val, alpha=0.4, color=c)
plt.ylim([25, 2100]) 
plt.xticks(rotation=45, ha='right')
plt.ylabel ('# Variant Calls')
plt.savefig("RNAvsWES.png",format='png',dpi=150,bbox_inches='tight')



