library("TCGAbiolinks")


#Get all individuals that have Solid Tissue Normal RNA-Seq
stn_rs_query <- GDCquery(
    project = "TCGA-LUSC", 
    data.category = "Sequencing Reads", 
    data.format = "bam",
    workflow.type = "STAR 2-Pass Genome",
    experimental.strategy = "RNA-Seq",
    sample.type = "Solid Tissue Normal"
)

stn_rs_indv=substr(getResults(stn_rs_query, cols = "cases"), 1, 12)


#Check for matching primary tumor 
# stn_rs_indv

# #Get all individuals that have Solid Tissue Normal RNA-Seq
# stn_rs_query <- GDCquery(
#     project = "TCGA-LUSC", 
#     data.category = "Sequencing Reads", 
#     data.format = "bam",
#     workflow.type = "STAR 2-Pass Genome",
#     experimental.strategy = "RNA-Seq",
#     sample.type = "Primary Tumor",
#     barcode = c(stn_indv)
# )

# stn_rs_indv=substr(getResults(stn_rs_query, cols = "cases"), 1, 12)

# stn_rs_indv






#Get those that also have somatic vcf made from Solid Tissue Normal vs primary tumor
somaticSNP_rs_query <- GDCquery(
    project = "TCGA-LUSC",
    data.category = "Simple Nucleotide Variation",
    data.format = "vcf",
    workflow.type = "MuTect2 Annotation",
    experimental.strategy = "WXS",
    barcode = c(stn_rs_indv),
    sample.type = "Solid Tissue Normal"
)

write.table(getResults(somaticSNP_rs_query),"somaticSNP.metadata.tsv" ,sep = '\t',row.names = FALSE,quote=FALSE)

somaticSNP_rs_indv=substr(getResults(somaticSNP_rs_query, cols = "cases"), 1, 12)





#Download RNA-Seq and VCF. Get rna-seq metadata. 
rs_download <- GDCquery(
    project = "TCGA-LUSC", 
    data.category = "Sequencing Reads", 
    data.format = "bam",
    workflow.type = "STAR 2-Pass Genome",
    experimental.strategy = "RNA-Seq",
    barcode = c(somaticSNP_rs_indv)
)

write.table(getResults(rs_download),"rs.metdata.tsv" ,sep = '\t',row.names = FALSE,quote=FALSE)

GDCdownload(
    query = rs_download, 
    method = "client",
    token.file="gdc-user-token.2022-11-14T21_18_50.458Z.txt",
    files.per.chunk=10
    )

GDCdownload(
    query = somaticSNP_rs_query, 
    method = "client",
    token.file="gdc-user-token.2022-11-14T21_18_50.458Z.txt",
    files.per.chunk=10
    )
