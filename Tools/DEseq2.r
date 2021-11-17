library(DESeq2)

#Read in count information. 
countData = read.table("JH-2021-6-11-Counts-Mason-Nasal.tsv",header=TRUE,row.names=NULL,sep = '\t')
#Remove duplicates
countData=countData[!duplicated(countData[ , c("Gene_stable_ID")]),]
##Read in expermental design
colData = read.table("Design",header=TRUE,row.names=NULL,sep = '\t')
colData=colData[colData$RT_PCR_Result=="Detected",]



#Get sample names 
samples=as.character(colData$samples)
#Extract samples from counts
counts=countData[,samples]
#Round to nearest int
counts=round(counts,0)
#Set Gene_stable_ID as rownames
rownames(counts)=countData$Gene_stable_ID




##Bring everything together 
dds = DESeqDataSetFromMatrix(countData = counts,colData = colData,design = ~ Ancestry)

dds <- estimateSizeFactors(dds)
idx <- rowSums( counts(dds, normalized=TRUE) >= 30 ) >= 10
dds <- dds[idx,]
dds <- DESeq(dds)


result = results(dds, contrast=c("Ancestry","EUR","AFR"))
result = result[complete.cases(result),]## to remove rows with NA
write.table(result, "EUR_AFR_DEG.txt")


result = results(dds, contrast=c("Ancestry","EUR","SAS"))
result = result[complete.cases(result),]## to remove rows with NA
write.table(result, "EUR_SAS_DEG.txt")

result = results(dds, contrast=c("Ancestry","EUR","EAS"))
result = result[complete.cases(result),]## to remove rows with NA
write.table(result, "EUR_EAS_DEG.txt")

result = results(dds, contrast=c("Ancestry","EUR","AMR"))
result = result[complete.cases(result),]## to remove rows with NA
write.table(result, "EUR_AMR_DEG.txt")






result = results(dds, contrast=c("Ancestry","EAS","AFR"))
result = result[complete.cases(result),]## to remove rows with NA
write.table(result, "EAS_AFR_DEG.txt")

result = results(dds, contrast=c("Ancestry","SAS","AFR"))
result = result[complete.cases(result),]## to remove rows with NA
write.table(result, "SAS_AFR_DEG.txt")

result = results(dds, contrast=c("Ancestry","AMR","AFR"))
result = result[complete.cases(result),]## to remove rows with NA
write.table(result, "AMR_AFR_DEG.txt")



result = results(dds, contrast=c("Ancestry","SAS","EAS"))
result = result[complete.cases(result),]## to remove rows with NA
write.table(result, "SAS_EAS_DEG.txt")

result = results(dds, contrast=c("Ancestry","SAS","AMR"))
result = result[complete.cases(result),]## to remove rows with NA
write.table(result, "SAS_AMR_DEG.txt")

result = results(dds, contrast=c("Ancestry","EAS","AMR"))
result = result[complete.cases(result),]## to remove rows with NA
write.table(result, "EAS_AMR_DEG.txt")
















