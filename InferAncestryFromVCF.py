import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn import metrics
import random
from pandas import DataFrame
import numpy as np
from os.path import exists
from sklearn.neural_network import MLPClassifier
from sklearn.preprocessing import StandardScaler
import os
import glob
# This code will suppress warnings.
import warnings
from sklearn.exceptions import ConvergenceWarning
warnings.simplefilter("ignore", ConvergenceWarning)

#Read in config file
configfile: "config.yaml"
Interval = config['Interval']
PCs=config['PCs']
OutputDir=config['OutputDir']

#ML methods##############

def neural_network(xtrain,ytrain,xtest,ytest,classification,unktest,sample_id,run):
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

    pred = models[idx].predict(unktest)
    pred=pd.DataFrame(data=pred.flatten())
    pred.columns =[classification]
    pred.insert(len(pred.columns),"TestAcc",[test_acc],True)
    pred.insert(len(pred.columns),"run_accession",[sample_id],True)



    if not os.path.isfile(classification+run+'NeuralNetworkResults'):
        pred.to_csv(classification+run+'NeuralNetworkResults',sep='\t',mode='a',index=False,header=True)
    else: # else it exists so append without writing the header
        pred.to_csv(classification+run+'NeuralNetworkResults',sep='\t',mode='a',index=False,header=None)

    pred.to_csv(OutputDir+ "/" + sample_id + "/" + classification+run+'NeuralNetworkResults',sep='\t',mode='w',index=False,header=True)

def RandomForest_Classifier(xtrain,ytrain,xtest,ytest,classification,unktest,sample_id,run):

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




    pred = models[idx].predict(unktest)
    pred=pd.DataFrame(data=pred.flatten())
    pred.columns =[classification]
    pred.insert(len(pred.columns),"TestAcc",[testacc],True)
    pred.insert(len(pred.columns),"run_accession",[sample_id],True)




    if not os.path.isfile(classification+run+'RandomForestResults'):
        pred.to_csv(classification+run+'RandomForestResults',sep='\t',mode='a',index=False,header=True)
    else: # else it exists so append without writing the header
        pred.to_csv(classification+run+'RandomForestResults',sep='\t',mode='a',index=False,header=None)

    pred.to_csv(OutputDir+ "/" + sample_id + "/" + classification+run+'RandomForestResults',sep='\t',mode='w',index=False,header=True)


#SVMs with rbf Kernels
def svm(xtrain,ytrain,xtest,ytest,classification,unktest,sample_id,run):
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
    pred.insert(len(pred.columns),"run_accession",[sample_id],True)


    df = pd.concat([pred, prob], axis=1)


    if not os.path.isfile(classification+run+'SVMResults'):
        df.to_csv(classification+run+'SVMResults',sep='\t',mode='a',index=False,header=True)
    else: # else it exists so append without writing the header
        df.to_csv(classification+run+'SVMResults',sep='\t',mode='a',index=False,header=None)

    df.to_csv(OutputDir+ "/" + sample_id + "/" + classification+run+'SVMResults',sep='\t',mode='w',index=False,header=True)


#N PCs
pcs=[]
for i in range(1,PCs+1):
    pcs.append("PC"+str(i))


def PCA(vcf,metadata_path,run,path,sample_id,classification):
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
        samples.loc[len(samples.index)] = [sample_id]



        #Specify all samples for PCA. This is needed for plink.
        file = pd.concat([samples, samples,samples], axis=1)
        file.to_csv(path+"samples", sep=" ",header=False,index=False)
        #Specify samples to train PCA with, all other samples will be projected.
        train=train.to_frame()
        train.columns =['Sample name']
        train.to_csv(path+"train", sep=" ",header=False,index=False)


        #shell("plink --vcf " + vcf + " --pca "+str(PCs)+" --within " + path + "samples --mac 1 --pca-clusters " + path + "train --out " + path + run)
        #Run PLINK PCA
        plink=open(path +"plink.sh", "w")
        plink.write("source activate Ancestry")
        plink.write('\n')
        plink.write("plink --vcf " + vcf + " --pca "+str(PCs)+" --within " + path + "samples --mac 1 --pca-clusters " + path + "train --out " + path + run)
        plink.close()
        subprocess.run("bash "+path +"plink.sh", shell=True, check=True)

        #PCA data from plink
        vardata= pd.read_csv(path + run + '.eigenvec',delimiter=r"\s+",names=['Drop','Sample name']+pcs)
        vardata=vardata.drop(['Drop'], axis=1)




        #Merge sample metadata with PCA data
        data=pd.merge(metadata,vardata,on=['Sample name'])



        #Make unkown, training, and test sets for ML methods.
        Train = pd.merge(data,train,on=['Sample name'])
        x_train=Train[pcs]  # Features
        y_train=Train[classification + ' code']  # Labels

        Test = pd.merge(data,test,on=['Sample name'])
        x_test=Test[pcs]  # Features
        y_test=Test[classification + ' code']  # Labels



        unk_test = vardata[vardata['Sample name'] == sample_id ]
        unk_test=unk_test[pcs]



        svm(x_train,y_train,x_test,y_test,classification,unk_test,sample_id,run)
        RandomForest_Classifier(x_train,y_train,x_test,y_test,classification,unk_test,sample_id,run)
        neural_network(x_train,y_train,x_test,y_test,classification,unk_test,sample_id,run)

#####################


#Create output directory



#Read in run accession IDs from txt file.
with open ("ids.txt") as f:
    id=f.read().splitlines()

chrPCSpec="Chr"+str(Interval).replace("[", "").replace("]", "").replace(",", "_").replace(" ", "")+".PC"+str(PCs)

if Interval == "All":
    chr_list = list(range(1, 23))
else:
    chr_list=Interval




rule all:
    input: expand("{wd}/{sample}/Superpopulation{chrPCSpec}SVMResults",sample=id,wd=OutputDir,chrPCSpec=chrPCSpec)


#File must be indexed
rule Intersect:
    input:
        "{wd}/{sample}/{sample}.vcf.gz"
    output:
        "{wd}/{sample}/Chr{chr}.{sample}.vcf.gz.tbi"

    shell:
        """
        bcftools view -r {wildcards.chr} --threads 7 --output-type z  {input} > {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.vcf.gz
        bcftools index -t {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.vcf.gz
        

        #Extract coordinates for input vcf  from 1KGP
        bcftools view -R {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.vcf.gz --threads 7 --output-type z data/Chr{wildcards.chr}_SNPs.vcf.gz > {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.1KGP.vcf.gz   
    
               
        #Remove headers and intersect varients for matching CHROM,POS,REF. ALT in gatk must be . or same as ALT in 1KGP. Outputs list of coordinates to use for filtering.
        gunzip -c    {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.1KGP.vcf.gz | grep -v "##" > {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.1KGP.vars
        gunzip -c   {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.GATK.vcf.gz | grep -v "##" > {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.GATK.vars
        """)
               
        gatk=pd.read_csv(wildcards.wd+'/'+wildcards.sample+'/Chr'+wildcards.chr+'.GATK.vars',sep='\t')                    
        KGP=pd.read_csv(wildcards.wd+'/'+wildcards.sample+'/Chr'+wildcards.chr+'.1KGP.vars',sep='\t')
        cords=gatk[(gatk['#CHROM'] == KGP['#CHROM']) & (gatk['POS'] == KGP['POS']) & (gatk['REF'] == KGP['REF']) &  ((gatk['ALT'] == KGP['ALT']) |  (gatk['ALT'] == '.' ) )][["#CHROM","POS"]] 
        cords.to_csv(wildcards.wd+'/'+wildcards.sample+'/Chr'+wildcards.chr+'Finalcords',sep='\t',index=False,header=None,mode='w')

        shell("""
        rm {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.1KGP.vars {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.GATK.vars
        
        #Filter
        bcftools view -R {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}Finalcords --threads 7 --output-type z {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.1KGP.vcf.gz > {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.1KGP.filtered.vcf.gz
        bcftools index -t {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.1KGP.filtered.vcf.gz
        bcftools view -R {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}Finalcords --threads 7 --output-type z {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.GATK.vcf.gz  > {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.GATK.filtered.vcf.gz
        bcftools index -t {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.GATK.filtered.vcf.gz
        rm {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.1KGP.vcf.gz {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.GATK.vcf.gz
        
        #Merge
        bcftools merge {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.1KGP.filtered.vcf.gz {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.GATK.filtered.vcf.gz --output-type z --threads 7 > {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.{wildcards.sample}.vcf.gz
        bcftools index -t {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.{wildcards.sample}.vcf.gz
        """)
 
    
 
        """

rule concat:
    input:
        expand("{wd}/{sample}/Chr{chr}.{sample}.vcf.gz.tbi",sample=id,chr=chr_list,wd=OutputDir)
    output:
        "{wd}/{sample}/"+chrPCSpec.replace(".PC"+str(PCs),"_")+".{sample}.vcf.gz"

    run:

        bcf_input=str(list(filter(lambda x: wildcards.sample in x, input))).replace("[", "").replace("]", "").replace(",", "").replace("'", "").replace(".tbi", "")
        path=wildcards.wd + "/"+wildcards.sample + "/"
        #shell("bcftools concat "+str(input).replace(".tbi", "")+" --output-type z --threads 7 > {wildcards.wd}/{wildcards.sample}/"+chrPCSpec.replace(".PC"+str(PCs),"_")+".{wildcards.sample}.vcf.gz")
        bcftools=open(path +"bcftools.sh", "w")
        bcftools.write("source activate Ancestry")
        bcftools.write('\n')
        bcftools.write("bcftools concat "+bcf_input+" --output-type z --threads 7 > "+ wildcards.wd+"/"+wildcards.sample+"/"+chrPCSpec.replace(".PC"+str(PCs),"_")+"."+wildcards.sample+".vcf.gz")
        bcftools.close()
        subprocess.run("bash "+path +"bcftools.sh", shell=True, check=True)
        
rule SuperPop:
    input: "{wd}/{sample}/"+chrPCSpec.replace(".PC"+str(PCs),"_")+".{sample}.vcf.gz"

    output:
        "{wd}/{sample}/Superpopulation{chrPCSpec}SVMResults"

    run:


        path=wildcards.wd + "/"+wildcards.sample + "/"
        vcf=str(input)
        metadata_path="data/1KGP.metadata.tsv"

        #All is ChrAll,Chr#,ect
        PCA(vcf,metadata_path,chrPCSpec,path,wildcards.sample,"Superpopulation")
