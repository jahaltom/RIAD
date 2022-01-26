import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn import metrics
import random
from pandas import DataFrame
import numpy as np
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import StandardScaler
import os
# This code will suppress warnings.
import warnings
from sklearn.exceptions import ConvergenceWarning
warnings.simplefilter("ignore", ConvergenceWarning)


#ML methods##############

def neural_network(xtrain,ytrain,xtest,ytest,classification,unktest,srrid,run):
    # Now make plots of Training Accuracy and plots of Testing Accuracy for NN, one layer 
    models=[]
    testacc = []
    bestlearnrate=[]
    best_alpha=[]
    noderange = range(1,11,1)
    
    alpharange = np.logspace(-6,0,5)
    learnrate = np.logspace(-3,-1,3)
    
    for nodecnt in noderange:
            
        bestacc = 0
        bestalpha = 0
        bestrate = 0   
        for alpha in alpharange:
            for rate in learnrate:
    
  
    
                clf = MLPClassifier(hidden_layer_sizes=(nodecnt),\
                            activation='relu', solver='sgd', \
                            learning_rate='adaptive', batch_size='auto',\
                            shuffle=True,alpha=alpha, \
                            nesterovs_momentum=True,verbose=False,\
                            momentum=0.9,early_stopping=False,\
                            max_iter=200, learning_rate_init=rate)
        
                clf.fit(xtrain, ytrain)
    
                if clf.score(xtest,ytest)> bestacc:
                    bestacc = clf.score(xtest,ytest)
                    bestalpha = alpha
                    bestrate = rate
                    mod=clf
    
      
     
    
        testacc.append(bestacc)
        bestlearnrate.append(bestrate)
        best_alpha.append(bestalpha)
        models.append(mod)
        
       
  
    
    idx = testacc.index(max(testacc))
  
   
    
    test_acc= models[idx].score(xtest,ytest)
    #Predict unk
    print(test_acc,max(testacc))
    #Predict unk
 
    pred = models[idx].predict(unktest)    
    pred=pd.DataFrame(data=pred.flatten())
    pred.columns =[classification]
    pred.insert(len(pred.columns),"TestAcc",[test_acc],True)
    pred.insert(len(pred.columns),"run_accession",[srrid],True)
    
   
        
    if not os.path.isfile(classification+run+'NeuralNetworkResults'):
        pred.to_csv(classification+run+'NeuralNetworkResults',sep='\t',mode='a',index=False,header=True)      
    else: # else it exists so append without writing the header
        pred.to_csv(classification+run+'NeuralNetworkResults',sep='\t',mode='a',index=False,header=None)    
    
    pred.to_csv(DIR+ "/" + srrid + "/" + classification+run+'NeuralNetworkResults',sep='\t',mode='w',index=False,header=True)    

def RandomForest_Classifier(xtrain,ytrain,xtest,ytest,classification,unktest,srrid,run):
     
    models=[]
    
    
    numtreerange = [1,5,10,25,50,100,200]
    overalltestacc = []
    overalltrainacc=[]
    bestNumTrees=[]
    depthrange = range(1,11)
    for depth in depthrange:
        
        best_num_trees = 0
        bestTest_acc = 0
        bestTrain_acc=0
       
        for num_trees in numtreerange:
            clf = RandomForestClassifier(bootstrap=True,n_estimators=num_trees,max_features='auto',criterion='gini',random_state=1,max_depth=depth)
            clf.fit(xtrain, ytrain)
           
            testacc= clf.score(xtest,ytest)
            trainacc=clf.score(xtrain,ytrain)
            if testacc >= bestTest_acc:
                bestTest_acc = testacc
                best_num_trees = num_trees
                bestTrain_acc=trainacc
                mod=clf



        models.append(mod)
        overalltrainacc.append( bestTrain_acc )
        overalltestacc.append(bestTest_acc) 
        bestNumTrees.append(best_num_trees)
        
    idx = overalltestacc.index(max(overalltestacc))
    #clf = RandomForestClassifier(bootstrap=True,n_estimators=bestNumTrees[idx],max_features=None,criterion='gini',max_depth=idx+1)
    #models[idx].fit(xtrain, ytrain)   
            
    testacc= models[idx].score(xtest,ytest)
    
    #Predict unk
    print(testacc,max(overalltestacc))
    
    
    pred = models[idx].predict(unktest)    
    pred=pd.DataFrame(data=pred.flatten())
    pred.columns =[classification]
    pred.insert(len(pred.columns),"TestAcc",[testacc],True)
    pred.insert(len(pred.columns),"run_accession",[srrid],True)
    
    
    
   
    if not os.path.isfile(classification+run+'RandomForestResults'):
        pred.to_csv(classification+run+'RandomForestResults',sep='\t',mode='a',index=False,header=True)    
    else: # else it exists so append without writing the header
        pred.to_csv(classification+run+'RandomForestResults',sep='\t',mode='a',index=False,header=None)    
    
    pred.to_csv(DIR+ "/" + srrid + "/" + classification+run+'RandomForestResults',sep='\t',mode='w',index=False,header=True)
   

#SVMs with Polynomial Kernels
def svm(xtrain,ytrain,xtest,ytest,classification,unktest,srrid,run):
    max_iter=1000    
    Cvals = np.logspace(-4,2,10)
    gamma_vals = np.logspace(-3, 2, 10)   
    

    bestc = 0
    bestgamma = 0
    besttrainacc = 0
    besttestacc = 0
    
    
   
    for c in Cvals:
        for gamma in gamma_vals:

            clf = SVC(C=c, kernel='rbf', gamma=gamma, random_state=1, shrinking=True, probability=True,max_iter=max_iter)

            clf.fit(xtrain, ytrain)

            trainacc = clf.score(xtrain, ytrain)
            testacc = clf.score(xtest, ytest)



            if testacc > besttestacc:
                besttestacc = testacc
                besttrainacc = trainacc
                bestgamma = gamma
                bestc = c
                mod=clf


    
    #Predict unk
        
    prob= mod.predict_proba(unktest)
    prob= pd.DataFrame(prob, columns = mod.classes_)

    
    pred= mod.predict(unktest)
    pred= pd.DataFrame(data=pred.flatten())
    pred.columns =[classification]
    
    pred.insert(len(pred.columns),"TestAcc",[besttestacc],True)
    pred.insert(len(pred.columns),"run_accession",[srrid],True)
    print(pred,prob)
    
    df = pd.concat([pred, prob], axis=1)
    
    
    if not os.path.isfile(classification+run+'SVMResults'):
        df.to_csv(classification+run+'SVMResults',sep='\t',mode='a',index=False,header=True) 
    else: # else it exists so append without writing the header
        df.to_csv(classification+run+'SVMResults',sep='\t',mode='a',index=False,header=None) 
    
    df.to_csv(DIR+ "/" + srrid + "/" + classification+run+'SVMResults',sep='\t',mode='w',index=False,header=True)


    
  
    
def PCA(vcf,metadata_path,run,path,srrid,classification):
        #Import metadata and gather sample names for each person in 1KGP.
        metadata=pd.read_csv(metadata_path,sep='\t')
        samples=metadata['Sample name'].copy()
        
        #Randomly suffle samples
        random.shuffle(samples)
        #Devide samples up into 0.70 train and 0.30 test. 
        x=len(samples)
        tr=x*0.70
        tr = int(round(tr, 0))
        tst=x-tr
        train=samples[0:tr]
        test=samples[tr:len(samples)]
        
        #Make test into dataframe and add header
        test=test.to_frame()
        test.columns =['Sample name']
      
        #Add unkown to 1KGP samples
        samples=samples.to_frame()
        samples.loc[len(samples.index)] = ['sample']
        
        
        
        #Specify all samples for PCA. This is needed for plink.  
        file = pd.concat([samples, samples,samples], axis=1)
        file.to_csv(path+"samples", sep=" ",header=False,index=False)
        #Specify samples to train PCA with, all other samples will be projected. 
        train=train.to_frame()
        train.columns =['Sample name']
        train.to_csv(path+"train", sep=" ",header=False,index=False)
        
        
        #Run PLINK PCA
        shell("plink --vcf " + vcf + " --pca 20 --within " + path + "samples --mac 1 --pca-clusters " + path + "train --out " + path + run)
        
        #PCA data from plink
        vardata= pd.read_csv(path + run + '.eigenvec',delimiter=r"\s+",names=['Drop','Sample name','PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10','PC11','PC12','PC13','PC14','PC15','PC16','PC17','PC18','PC19','PC20'])
        vardata=vardata.drop(['Drop'], axis=1)
        
        
  
        
        #Merge sample metadata with PCA data
        data=pd.merge(metadata,vardata,on=['Sample name'])
        

        
        #Make unkown, training, and test sets for ML methods.
        Train = pd.merge(data,train,on=['Sample name'])
        x_train=Train[['PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10','PC11','PC12','PC13','PC14','PC15','PC16','PC17','PC18','PC19','PC20']]  # Features
        y_train=Train[classification + ' code']  # Labels
        
        Test = pd.merge(data,test,on=['Sample name'])
        x_test=Test[['PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10','PC11','PC12','PC13','PC14','PC15','PC16','PC17','PC18','PC19','PC20']]  # Features
        y_test=Test[classification + ' code']  # Labels
        
        
        
        unk_test = vardata[vardata['Sample name'] == "sample" ]
        unk_test=unk_test[['PC1','PC2','PC3','PC4','PC5','PC6','PC7','PC8','PC9','PC10','PC11','PC12','PC13','PC14','PC15','PC16','PC17','PC18','PC19','PC20']]
      
    
            
        svm(x_train,y_train,x_test,y_test,classification,unk_test,srrid,run)
        #RandomForest_Classifier(x_train,y_train,x_test,y_test,classification,unk_test,srrid,run)
        #neural_network(x_train,y_train,x_test,y_test,classification,unk_test,srrid,run)

      


#Read in run accession IDs from txt file. 
with open ("RAids.txt") as f:
    ra=f.read().splitlines()
# output directory
DIR='output'

test="ChrAll_PC20"




phyla="FILE"

rule all:
	input: ["{wd}/{sample}/Population{test}SVMResults".format(sample=s,wd=DIR,test=test) for s in ra]

rule all:
    input: expand("{wd}/{sample}/Superpopulation{test}SVMResults",sample=ra,wd=DIR,test=test)



        
rule SuperPop:
    input: "{wd}/{sample}/All.{sample}.vcf.gz"
    
    output:
        "{wd}/{sample}/Superpopulation{test}SVMResults"
    run:            
        
        #All.{sample}.vcf.gz must be indexted before bcftools can be used on it.
        shell("bcftools index -f -t {input}")
        
        
        #This first
        shell("tabix -h -R data/" + phyla +" {input} | bgzip  > {wildcards.wd}/{wildcards.sample}/" + phyla.replace(".tsv", ".vcf.gz"))

       
        srrid=str(output).split("/")[1]
        path=DIR + "/"+srrid + "/"
        vcf=path + phyla.replace(".tsv", ".vcf.gz")
 
        
        metadata_path="data/1KGP.metadata.tsv"    
        
        #All is ChrAll,Chr#,ect
        PCA(vcf,metadata_path,test,path,srrid,"Superpopulation")
     

        
        


