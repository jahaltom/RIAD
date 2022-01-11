import pandas as pd
from pandas import DataFrame


#Make list 1-22
chr_num=[]
for x in range(1,23):
      chr_num.append(str(x))
      
      
      


#Read in metadata
metadata=pd.read_csv("metadata",sep="\t")
#Gather BioProj_Population 
ids=metadata[['BioProj_Population']].drop_duplicates()
#Make list
ids_list=ids['BioProj_Population'].tolist() 
#df to store accuracy results for each Chr
results_Chr = pd.DataFrame()

for c in chr_num:
 
    #Read in ancestry inference results 
    results=pd.read_csv("SuperpopulationChr"+c+".PC20SVMResults",sep="\t")       
    
    #Merge metadata with results  
    resultsMeta=pd.merge(results,metadata,on=['run_accession'])               
    
    acclist=[]
    sample_size=[]
    numCorrect=0
    for i in ids_list:
       #Get sample size
       sampleSize=len(resultsMeta[resultsMeta['BioProj_Population'] == i])
       sample_size.append(sampleSize)
       #Calculate accuracy as % 
       correct=len(resultsMeta[(resultsMeta['BioProj_Population'] == i) & ((resultsMeta['Superpopulation'] == resultsMeta['Eth1']) | (resultsMeta['Superpopulation'] == resultsMeta['Eth2']))])
       acclist.append(correct/ sampleSize*100)   
       
       #Keep track of total corretly infered. 
       numCorrect=numCorrect + correct
       
    #Overall accuracy across studies
    acclist.append(numCorrect/sum(sample_size)*100)
    
    #Accuracy results for a Chr
    temp_results_Chr=DataFrame(acclist,columns=['Chr'+c])
    
    #Add to df that will store accuracy results for each Chr
    results_Chr = pd.concat([results_Chr, temp_results_Chr], axis=1)
    

#Add total sample size to sample_size and make into df
sample_size.append(sum(sample_size))
sample_size=DataFrame(sample_size,columns=['Sample Size'])  
#Add OverallAccuracy string to ids
ids.loc[-1] = ['OverallAccuracy']
ids=ids.reset_index(drop=True)
#Combine results
results_Chr=pd.concat([ids,sample_size,results_Chr], axis=1)


results_Chr.to_csv('Chr1-22_PC20_DPLT5.tsv',sep='\t',mode='w',index=False) 
