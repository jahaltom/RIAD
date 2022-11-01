import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sb
import os

os.chdir(r"C:\Users\15154\Documents")


#Read in results
results=pd.read_csv("AccuracyTable.tsv",index_col=0,sep="\t")

#List to store dfs
dfs = []
#List of columns wanting extraction
col_list=['All']
for c in range(1,23): 
    col_list.append(c)
    
    
#Loop through chrs 
for i in col_list:   
    #Extract single chr column 
    chrCol=[col for col in results.columns if 'Chr'+str(i)+".PC20SVMResults" in col]
    #Extract accuracies 
    df=results[chrCol].tail(7)
    #Add column for Chr num 
    df['Chr']=i
    #Rownames to column
    df.index.name = 'Ethnicity'
    df.reset_index(inplace=True)
    #Remove header
    df.columns = range(df.shape[1])
    #Drop Malaysian study
    df=df.drop([3, 3])
  
    dfs.append(df)

dfs = pd.concat(dfs)
dfs.columns =['Ethnicity', 'Accuracy', 'Chromosome']


plot=sb.barplot(x='Chromosome', y='Accuracy', hue="Ethnicity", data=dfs, ci = None)
plot.set_xticklabels(plot.get_xticklabels(), rotation=90,size = 7)  



plt.legend(bbox_to_anchor=(1.01, 1),borderaxespad=0)
plt.xlabel("Chromosome")
plt.ylabel("Accuracy")

#Keeps stuff being cropped


for patch in plot.patches :
    current_width = patch.get_width()
    diff = current_width - .1

    # we change the bar width
    patch.set_width(.1)

    # we recenter the bar
    patch.set_x(patch.get_x() + diff *.5)




plt.savefig("AccuracyByChr.png",format='png',dpi=150,bbox_inches='tight')
