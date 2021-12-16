import pandas as pd
from pandas import DataFrame



Super=['SuperpopulationChrAll.PC20SVMResults',  'SuperpopulationChrAll.PC20RandomForestResults',  'SuperpopulationChrAll.PC20NeuralNetworkResults']



for s in Super:
    #Read in metadata
    metadata=pd.read_csv("metadata",sep="\t")
    #Gather population names
    ids=metadata[['Population']].drop_duplicates()
    ids=ids['Population'].tolist()
    
    #Read in ancestry inference results 
    results=pd.read_csv(s,sep="\t")
    results=results.rename(columns = {'ID':'run_accession'})
        
    #Merge metadata with results  
    test=pd.merge(results,metadata,on=['run_accession'])
        
        
        
    acclist=[]
    sample_size=[]
    for i in ids:
       #Calculate accuracy as % 
       acclist.append(len(test[(test['Population'] == i) & ((test['Superpopulation'] == test['Eth1']) | (test['Superpopulation'] == test['Eth2']))]) / len(test[test['Population'] == i])*100)
       #Get whole sample size
       sample_size.append(len(test[test['Population'] == i]))
        
        
        
    #Gather population names as df
    ids=metadata[['Population']].drop_duplicates()
    ids=ids.reset_index(drop=True)
    
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
