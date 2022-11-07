import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sb


os.chdir(r"C:\Users\15154\Documents")


################
Ancestry_Results_File='SuperpopulationChrAll.PC20SVMResults'
outfilePrefix="RIA"

#################

#Read in metadata
metadata=pd.read_csv("metadata.txt",sep="\t")

#Read in ancestry inference results 
results=pd.read_csv(Ancestry_Results_File,sep="\t")       

#Merge metadata with results  
resultsMeta=pd.merge(results,metadata,on=['run_accession'])  

#Set whole field to No then change based off accuracy 
resultsMeta['Accurate'] ="No"
for idn in resultsMeta.index:
    if ((resultsMeta['Superpopulation'][idn] == resultsMeta['Eth1'][idn]) | (resultsMeta['Superpopulation'][idn] == resultsMeta['Eth2'][idn]))==True:
        resultsMeta['Accurate'][idn] ="Yes"

#Group self-reported ethnicity and accuracy 
resultsMeta["Accurate_Eth"]=resultsMeta["Accurate"]+resultsMeta["Eth1"]

##Num of SNPs
#Add in dummy row to sperate accurate from non.  
resultsMeta = pd.concat([resultsMeta[resultsMeta.Accurate == 'No'], pd.DataFrame({'Num of SNPs': [0,0,0], 'Accurate_Eth': [""," ","  "], 'Eth1':["","",""]}), resultsMeta[resultsMeta.Accurate == 'Yes']])
#Plot
plot=sb.scatterplot(resultsMeta["Num of SNPs"],resultsMeta["Accurate_Eth"], hue=resultsMeta['Eth1'], data=resultsMeta,hue_order=["","AFR","AMR","EAS","EUR","SAS"],palette=["white","red","green","blue","purple","orange"])
plt.margins(y=0.1)
plt.xlabel("Total SNPs/Sample")
plt.ylabel("Ancestry Correctly Inferred")

plt.legend(title='Superpopulation')
plt.grid(True, axis='x')
plt.savefig(outfilePrefix+"AccuracyBySNP#.png",format='png',dpi=150,bbox_inches='tight')


plt.show()


##Num of Unique Mapped Reads
#Add in dummy row to sperate accurate from non.  
resultsMeta = pd.concat([resultsMeta[resultsMeta.Accurate == 'No'], pd.DataFrame({'Num of Unique Mapped Reads': [0,0,0], 'Accurate_Eth': [""," ","  "], 'Eth1':["","",""]}), resultsMeta[resultsMeta.Accurate == 'Yes']])
#Plot
plot=sb.scatterplot(resultsMeta["Num of Unique Mapped Reads"],resultsMeta["Accurate_Eth"], hue=resultsMeta['Eth1'], data=resultsMeta,hue_order=["","AFR","AMR","EAS","EUR","SAS"],palette=["white","blue","orange","green","red","purple"])
plt.margins(y=0.1)
plt.xlabel("# of Unique Mapped Reads/Sample")
plt.ylabel("Ancestry Correctly Inferred")

plt.legend(title='Superpopulation')
plt.grid(True, axis='x')
plt.savefig(outfilePrefix+"AccuracyBy#UniqMappedReads.png",format='png',dpi=150,bbox_inches='tight')


