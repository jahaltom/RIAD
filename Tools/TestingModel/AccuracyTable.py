import pandas as pd
from pandas import DataFrame

Super=['SuperpopulationChrAll.PC20SVMResults',  'SuperpopulationChrAll.PC20RandomForestResults',  'SuperpopulationChrAll.PC20NeuralNetworkResults']


#Read in metadata
metadata=pd.read_csv("metadata",sep="\t")
#Gather BioProj_Population 
ids=metadata[['BioProj_Population']].drop_duplicates()


#Make list
ids_list=ids['BioProj_Population'].tolist()  

for s in Super:
  
    
    #Read in ancestry inference results 
    results=pd.read_csv(s,sep="\t")       
    
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
    
   
    
    if s.endswith("NeuralNetworkResults"):
        resultsNN=DataFrame(acclist,columns=['Neural Network Accuracy'])
    elif s.endswith("RandomForestResults"):
        resultsRF=DataFrame(acclist,columns=['Random Forest Accuracy'])      
    else:
        resultsSVM=DataFrame(acclist,columns=['SVM Accuracy'])

#Add total sample size to sample_size and make into df
sample_size.append(sum(sample_size))
sample_size=DataFrame(sample_size,columns=['Sample Size'])  
#Add OverallAccuracy string to ids
ids.loc[-1] = ['OverallAccuracy']
ids=ids.reset_index(drop=True)
#Combine results
accResults=pd.concat([ids,sample_size,resultsRF,resultsNN,resultsSVM], axis=1)


accResults.to_csv('ChrAll_PC15_DPLT5.tsv',sep='\t',mode='w',index=False) 
