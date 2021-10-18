library(DESeq2)

#Read in count information. 
countData = read.table("JH-2021-6-11-Counts-Mason-Nasal.tsv",header=TRUE,row.names=NULL,sep = '\t')
#Remove duplicates
countData=countData[!duplicated(countData[ , c("Gene_stable_ID")]),]
##Read in expermental design
colData = read.table("Design",header=TRUE,row.names=NULL,sep = '\t')

#Get sample names 
samples=as.character(colData$samples)
#Extract samples from counts
counts=countData[,samples]
#Round to nearest int
counts=round(counts,0)
#Set Gene_stable_ID as rownames
rownames(counts)=countData$Gene_stable_ID



##Bring everything together 
dds = DESeqDataSetFromMatrix(countData = counts,colData = colData,design = ~ Ancestry + RT_PCR_Result)

dds <- estimateSizeFactors(dds)
idx <- rowSums( counts(dds, normalized=TRUE) >= 5 ) >= 3
dds <- dds[idx,]
dds <- DESeq(dds)


result - results(dds, contrast=c(Ancestry,SAS,EAS))
result - result[complete.cases(result),]## to remove rows with NA
write.table(result, TumorDEG.txt)
