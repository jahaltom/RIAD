import pandas as pd
from pandas import DataFrame



Super=['SuperpopulationChrAll.PC20SVMResults',  'SuperpopulationChrAll.PC20RandomForestResults',  'SuperpopulationChrAll.PC20NeuralNetworkResults']

#Read in metadata
metadata=pd.read_csv("metadata",sep="\t")

for s in Super:
 
    #Gather BioProj_Population 
    ids=metadata[['BioProj_Population']].drop_duplicates()
    ids=ids.reset_index(drop=True)
    #Make list
    ids_list=ids['BioProj_Population'].tolist()
    
    #Read in ancestry inference results 
    results=pd.read_csv(s,sep="\t")
    results=results.rename(columns = {'ID':'run_accession'})
        
    #Merge metadata with results  
    resultsMeta=pd.merge(results,metadata,on=['run_accession'])
        
        
       
    acclist=[]
    sample_size=[]
    for i in ids_list:
       #Calculate accuracy as % 
       acclist.append(len(resultsMeta[(resultsMeta['BioProj_Population'] == i) & ((resultsMeta['Superpopulation'] == resultsMeta['Eth1']) | (resultsMeta['Superpopulation'] == resultsMeta['Eth2']))]) / len(resultsMeta[resultsMeta['BioProj_Population'] == i])*100)
       #Get sample size
       sample_size.append(len(resultsMeta[resultsMeta['BioProj_Population'] == i]))
        
        
        
    
    #Make sample_size into df
    sample_size=DataFrame(sample_size,columns=['Sample Size'])
    
    
    if s.endswith("NeuralNetworkResults"):
        resultsNN=DataFrame (acclist,columns=['Neural Network Accuracy'])
    elif s.endswith("RandomForestResults"):
        resultsRF=DataFrame (acclist,columns=['Random Forest Accuracy'])      
    else:
        resultsSVM=DataFrame (acclist,columns=['SVM Accuracy'])
    
    





accResults=pd.concat([ids,sample_size,resultsRF,resultsNN,resultsSVM], axis=1)
accResults.to_csv('ChrAll_PC20_DPLT5',sep='\t',mode='w',index=False) 
