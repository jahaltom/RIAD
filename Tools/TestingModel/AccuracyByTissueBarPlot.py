import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sb
import os

os.chdir(r"C:\Users\15154\Documents")

#Read in metadata
metadata=pd.read_csv("metadata.txt",sep="\t")
#Read in results
results=pd.read_csv("SuperpopulationChrAll.PC20SVMResults",sep="\t")



#Merge on run_accession
merged=pd.merge(metadata,results,on=['run_accession'])

#This column will be used to get sample size later
merged['sample size']=1

#adds 1 if ancestry inference accurate, 0 otherwise 
def test(s):
    if ((s['Superpopulation'] == s['Eth1']) | (s['Superpopulation'] == s['Eth2'])):
        return 1
    else:
        return 0
                                         
merged['Accuracy'] = merged.apply(test, axis=1)    
                                 
merged['TissueGroup']=merged['BioProj_Population'] + "#" + merged['Tissue']+"#"+merged['Eth1']+"#"+merged['Eth2']
                     

#Sum sample size and Accuracy
merged = merged.groupby(['TissueGroup'],as_index = False).sum()
#Calulate overall accuracy
merged['OverallAccuracy']=merged['Accuracy']/merged['sample size']*100



#Split TissueGroup column into 2
merged[['Study','Tissue','Eth1','Eth2']] = merged.TissueGroup.str.split("#",expand=True)


#Sort by self-reported ethnicity and study. 
merged=merged.sort_values(by=['Eth1',"Study"])



plot=sb.barplot(x='Tissue', y='OverallAccuracy', hue="Eth1", data=merged, ci = None)
plot.set_xticklabels(plot.get_xticklabels(), rotation=90,size = 7)  



plt.legend(bbox_to_anchor=(1.01, 1),borderaxespad=0)
plt.xlabel("Tissue")
plt.ylabel("Accuracy %")

#Keeps stuff being cropped


for patch in plot.patches :
    current_width = patch.get_width()
    diff = current_width - .1

    # we change the bar width
    patch.set_width(.2)

    # we recenter the bar
    patch.set_x(patch.get_x() + diff *.5)





plt.savefig("AccuracyByTissueEth.png",format='png',dpi=150,bbox_inches='tight')
