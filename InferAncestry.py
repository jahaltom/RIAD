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
        RandomForest_Classifier(x_train,y_train,x_test,y_test,classification,unk_test,srrid,run)
        neural_network(x_train,y_train,x_test,y_test,classification,unk_test,srrid,run)
        
#####################

        
#Create output directory
DIR='output'


#Read in run accession IDs from txt file.
with open ("RAids.txt") as f:
    ra=f.read().splitlines()

test="ChrAll.PC20"

chr_list = list(range(1, 23))


rule all:
    input: expand("{wd}/{sample}/Superpopulation{test}SVMResults",sample=ra,wd=DIR,test=test)


rule var_call:
    input:
        "{wd}/{sample}/2ndPass.Aligned.sortedByCoord.out.bam"
    output:
        "{wd}/{sample}/Chr{chr}.final.vcf.gz.tbi"
    group:
        "Genotype"
    shell:
        """       
        samtools view -b {input} {wildcards.chr}  > {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.bam



        gatk MarkDuplicates \
            I= {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.bam \
            O= {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.MarkDup.bam \
            CREATE_INDEX= true \
            METRICS_FILE= {wildcards.wd}/{wildcards.sample}/marked_dup_metrics.{wildcards.chr}.txt \
            VALIDATION_STRINGENCY= SILENT
        rm {wildcards.wd}/{wildcards.sample}/*Chr{wildcards.chr}.bam

        gatk AddOrReplaceReadGroups \
            I= {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.MarkDup.bam \
            O= {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.Grouped.bam \
            RGSM= sample \
            RGLB= lib \
            RGPL= plat \
            RGPU= plat
        rm {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.MarkDup*

        gatk SplitNCigarReads \
            -I {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.Grouped.bam  \
            -O {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.final.bam \
            -R data/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna
        rm {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.Grouped.bam



        gatk --java-options "-Xmx20g" HaplotypeCaller \
           -I {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.final.bam \
           -O {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.vcf.gz \
           -R data/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna \
           -L data/Chr{wildcards.chr}_SNPs.bed \
           -ERC GVCF \
           --native-pair-hmm-threads 7
        rm {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.final*


        gunzip -c {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.vcf.gz |  grep -v "0/0:0:0:0:0,0,0" | bgzip > {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.final.vcf.gz
        rm {wildcards.wd}/{wildcards.sample}/*Chr{wildcards.chr}.vcf.gz

        bcftools index -t {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.final.vcf.gz
        """
rule Genotype:
    input: 
        "{wd}/{sample}/Chr{chr}.final.vcf.gz.tbi"
    output: 
        "{wd}/{sample}/Chr{chr}.final2.vcf.gz.tbi"
    group: 
        "Genotype"
    shell:
        """
        gatk --java-options "-Xmx10g" GenotypeGVCFs \
        -R data/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna \
        -V {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.final.vcf.gz \
        -O {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.vcf.gz \
        --include-non-variant-sites
        

        gatk --java-options "-Xmx10g" VariantFiltration \
        -R data/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna \
        -V {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.vcf.gz \
        -O {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.filtered.vcf.gz \
        --filter-expression "DP < 5.0" \
        --filter-name "filter_DP" \
        --filter-expression "MQ < 40.0" \
        --filter-name "filter_MQ" \
        --filter-expression "FS > 60.0" \
        --filter-name "filter_FS" \
        --filter-expression "MQRankSum < -12.5" \
        --filter-name "filter_MQRankSum" \
        --filter-expression "ReadPosRankSum < -8.0" \
        --filter-name "filter_ReadPosRankSum" \
        --filter-expression "QD < 2.0" \
        --filter-name "filter_QD"
        rm {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.vcf.gz*


        gatk --java-options "-Xmx10g" SelectVariants \
        -R data/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna \
        -V {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.filtered.vcf.gz  \
        -O {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.final2.vcf.gz \
        --exclude-filtered true \
        --select-type-to-exclude INDEL

        rm {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.filtered.vcf.gz*
        """
rule compare:
    input: 
        "{wd}/{sample}/Chr{chr}.final2.vcf.gz.tbi"
    output: 
        "{wd}/{sample}/Chr{chr}.{sample}.vcf.gz.tbi"
    group:               
        "Genotype"
    shell:
        """

        gunzip -c  {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.final2.vcf.gz | grep -v  "\./\." | bgzip > {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.final2.filtered.vcf.gz
        bcftools index -t {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.final2.filtered.vcf.gz

        

        #all sites with matching positions. Will result in all 0/0 calls that match and will have 0/1 and 1/1 which also match but some will be the wronf alt allele.

        tabix -h -T data/Chr{wildcards.chr}_SNPs.bed {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.final2.filtered.vcf.gz | bgzip  > {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.comm_pos.vcf.gz
        bcftools index -t {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.comm_pos.vcf.gz


        rm {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.final2.filtered.vcf.gz*

        #Make an index file of comm_pos.vcf.gz. Handels blank vcfs
        gunzip -c {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.comm_pos.vcf.gz | {{ grep -v "#" || true; }} | awk -F'\t' -v OFS='\t' '{{print $1,$2}}' | cat > {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.comm_pos.interval_list

        #Extract positions from 1KGP
        tabix -h -R {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.comm_pos.interval_list  data/Chr{wildcards.chr}_SNPs.vcf.gz | bgzip   > {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}1000.vcf.gz
        bcftools index -t {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}1000.vcf.gz

        #Match only on positions that have the same allele.
        bcftools isec -c none {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.comm_pos.vcf.gz {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}1000.vcf.gz --output-type z --threads 7 -p {wildcards.wd}/{wildcards.sample}/dir{wildcards.chr}
        #0000.vcf.gz is records private to  comm_pos.vcf.gz.They will be 0/0 calls that are true matches and 1/1 and 0/1 calls with wrong alt allele. Remove 0/0 and same others as an index list.
        gunzip -c {wildcards.wd}/{wildcards.sample}/dir{wildcards.chr}/0000.vcf.gz | grep -v "0/0" | {{ grep -v "#" || true; }} | awk -F'\t' -v OFS='\t' '{{print $1,$2}}' | cat > {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}List

        #If List is not empty, then remove unwanted SNPs, else remove nothing.

        if test -s {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}List ; then
                bcftools view -T ^{wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}List {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}1000.vcf.gz --output-type z --threads 7 > {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.1KG.vcf.gz
                bcftools view -T ^{wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}List {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.comm_pos.vcf.gz --output-type z --threads 7 > {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.unk.vcf.gz
        else
                mv {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}1000.vcf.gz {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.1KG.vcf.gz
                mv {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.comm_pos.vcf.gz {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.unk.vcf.gz
        fi

        rm -r  {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.comm_pos.interval_list  {wildcards.wd}/{wildcards.sample}/dir{wildcards.chr} {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}1000.vcf.gz* {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.comm_pos.vcf.gz* {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}List

        #Use List to remove unmatched positions from data comtaining all matches 0/0 1/0 1/1.

        bcftools index -t {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.unk.vcf.gz
        bcftools index -t {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.1KG.vcf.gz
        #Combine 1KGP with sample to have ancestry infered

        bcftools merge {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.1KG.vcf.gz {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.unk.vcf.gz  --output-type z --threads 7 > {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.{wildcards.sample}.vcf.gz
        bcftools index -t {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.{wildcards.sample}.vcf.gz

        rm {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.1KG.vcf.gz* {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.unk.vcf.gz*

        """

rule concat:
    input: 
        expand("{wd}/{sample}/Chr{chr}.{sample}.vcf.gz.tbi",sample=ra,chr=chr_list,wd=DIR)
    output: 
        "{wd}/{sample}/All.{sample}.vcf.gz"
    group:
        "Infer Ancestry"
    shell:
        """
        rm {wildcards.wd}/{wildcards.sample}/2ndPass.Aligned.sortedByCoord.out.bam*
        rm {wildcards.wd}/{wildcards.sample}/*.final.vcf.gz*
        rm {wildcards.wd}/{wildcards.sample}/*.final2.vcf.gz*
        bcftools concat {wildcards.wd}/{wildcards.sample}/*{wildcards.sample}.vcf.gz --output-type z --threads 7 > {wildcards.wd}/{wildcards.sample}/All.{wildcards.sample}.vcf.gz
        """

rule SuperPop:
    input: "{wd}/{sample}/All.{sample}.vcf.gz"
    
    output:
        "{wd}/{sample}/Superpopulation{test}SVMResults"
    group:
        "Infer Ancestry"
    run:
   
        srrid=str(output).split("/")[1]
        path=DIR + "/"+srrid + "/"
        vcf=str(input)
        metadata_path="data/1KGP.metadata.tsv"    
        
        #All is ChrAll,Chr#,ect
        PCA(vcf,metadata_path,test,path,srrid,"Superpopulation")

        
