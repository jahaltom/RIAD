import pandas as pd
from pandas import DataFrame

Super=['SuperpopulationChrAll_PC20.Ann_Amniota.tsvSVMResults',
'SuperpopulationChrAll_PC20.Ann_Boreoeutheria.tsvSVMResults',
'SuperpopulationChrAll_PC20.Ann_Catarrhini.tsvSVMResults',
'SuperpopulationChrAll_PC20.Ann_Euarchontoglires.tsvSVMResults',
'SuperpopulationChrAll_PC20.Ann_Eutheria.tsvSVMResults',
'SuperpopulationChrAll_PC20.Ann_Haplorrhini.tsvSVMResults',
'SuperpopulationChrAll_PC20.Ann_Hominidae.tsvSVMResults',
'SuperpopulationChrAll_PC20.Ann_Homininae.tsvSVMResults',
'SuperpopulationChrAll_PC20.Ann_Homo_sapiens.tsvSVMResults',
'SuperpopulationChrAll_PC20.Ann_Primates.tsvSVMResults',
'SuperpopulationChrAll_PC20.Ann_Simiiformes.tsvSVMResults',
'SuperpopulationChrAll_PC20.EB_Amniota.tsvSVMResults',
'SuperpopulationChrAll_PC20.EB_Ann_Amniota.tsvSVMResults',
'SuperpopulationChrAll_PC20.EB_Ann_Boreoeutheria.tsvSVMResults',
'SuperpopulationChrAll_PC20.EB_Ann_Catarrhini.tsvSVMResults',
'SuperpopulationChrAll_PC20.EB_Ann_Euarchontoglires.tsvSVMResults',
'SuperpopulationChrAll_PC20.EB_Ann_Eutheria.tsvSVMResults',
'SuperpopulationChrAll_PC20.EB_Ann_Haplorrhini.tsvSVMResults',
'SuperpopulationChrAll_PC20.EB_Ann_Hominidae.tsvSVMResults',
'SuperpopulationChrAll_PC20.EB_Ann_Homininae.tsvSVMResults',
'SuperpopulationChrAll_PC20.EB_Ann_Homo_sapiens.tsvSVMResults',
'SuperpopulationChrAll_PC20.EB_Ann_Primates.tsvSVMResults',
'SuperpopulationChrAll_PC20.EB_Ann_Simiiformes.tsvSVMResults',
'SuperpopulationChrAll_PC20.EB_Boreoeutheria.tsvSVMResults',
'SuperpopulationChrAll_PC20.EB_Catarrhini.tsvSVMResults',
'SuperpopulationChrAll_PC20.EB_Euarchontoglires.tsvSVMResults',
'SuperpopulationChrAll_PC20.EB_Eutheria.tsvSVMResults',
'SuperpopulationChrAll_PC20.EB_Haplorrhini.tsvSVMResults',
'SuperpopulationChrAll_PC20.EB_Hominidae.tsvSVMResults',
'SuperpopulationChrAll_PC20.EB_Homininae.tsvSVMResults',
'SuperpopulationChrAll_PC20.EB_Homo_sapiens.tsvSVMResults',
'SuperpopulationChrAll_PC20.EB_Primates.tsvSVMResults',
'SuperpopulationChrAll_PC20.EB_Simiiformes.tsvSVMResults']


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
       if sampleSize > 0:
           #Calculate accuracy as % 
           correct=len(resultsMeta[(resultsMeta['BioProj_Population'] == i) & ((resultsMeta['Superpopulation'] == resultsMeta['Eth1']) | (resultsMeta['Superpopulation'] == resultsMeta['Eth2']))])
           acclist.append(correct/ sampleSize*100)   
           
           #Keep track of total corretly infered. 
           numCorrect=numCorrect + correct
       else:
           acclist.append("NA")
    #Overall accuracy across studies
    acclist.append(numCorrect/sum(sample_size)*100)
    
       
    ids_temp=ids.copy()
    
    resultsSVM=DataFrame(acclist,columns=['SVM Accuracy'])
    
    #Add total sample size to sample_size and make into df
    sample_size.append(sum(sample_size))
    sample_size=DataFrame(sample_size,columns=['Sample Size'])  
    #Add OverallAccuracy string to ids
    ids_temp.loc[-1] = ['OverallAccuracy']
    ids_temp=ids_temp.reset_index(drop=True)
    #Combine results
    accResults=pd.concat([ids_temp,sample_size,resultsSVM], axis=1)  
    accResults.to_csv(s+"DONE.tsv",sep='\t',mode='w',index=False) 
