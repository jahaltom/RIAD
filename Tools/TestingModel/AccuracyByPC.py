import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sb
import os

os.chdir(r"C:\Users\15154\Documents")


#Read in results
results=pd.read_csv("AccuracyTable.tsv",index_col=0,sep="\t")

#List to store dfs
dfs = []
#List of columns wanting extraction. PCs
col_list=[]
for pc in (2,5,10,15,20): 
    col_list.append(pc)
#Machine learning methods
ml_methods=["SVM","RandomForest","NeuralNetwork"]
#Loop through chrs 
for i in col_list: 
    for ml in ml_methods:
        #Extract single chr column 
        pcCol=[col for col in results.columns if 'ChrAll.PC'+str(i)+ml in col and 'Mask' not in col]
        #Extract accuracies 
        df=results[pcCol].loc[["OverallAccuracy"]]
    
        #Add column for PC and ML tool 
        df['PC']=i
        df['Machine Learning Method']=ml
        
        df.reset_index(inplace=True)
        
        #Remove header
        df.columns = range(df.shape[1]) 
        df=df.drop([0],axis=1)
        dfs.append(df)

    

dfs = pd.concat(dfs)
dfs.columns =['Accuracy', 'PC','Machine Learning Method']
print(df)

plot=sb.barplot(x='PC', y='Accuracy', hue="Machine Learning Method", data=dfs, ci = None)
plot.set_xticklabels(plot.get_xticklabels(), rotation=90,size = 7)  



plt.legend(bbox_to_anchor=(1.01, 1),borderaxespad=0)
plt.xlabel("Principal Components")
plt.ylabel("Accuracy %")

#Keeps stuff being cropped


for patch in plot.patches :
    current_width = patch.get_width()
    diff = current_width - .1

    # we change the bar width
    patch.set_width(.1)

    # we recenter the bar
    patch.set_x(patch.get_x() + diff *.5)




plt.savefig("AccuracyByPC.png",format='png',dpi=150,bbox_inches='tight')
