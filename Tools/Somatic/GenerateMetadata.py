import pandas as pd



#RNA-Seq TCGA metadata
rs_data=pd.read_csv("rs.metdata.tsv",sep='\t')
#Split and merge by individual
pt_rs_data=rs_data.loc[rs_data['sample_type'] == "Primary Tumor"]
stn_rs_data=rs_data.loc[rs_data['sample_type'] == "Solid Tissue Normal"]
rs_data=pd.merge(pt_rs_data,stn_rs_data,on=["cases.submitter_id"])
#id_x and id_y contain the names of the directories that contain a single RNA-Seq sample for tumor and normal. 
rs_data = rs_data.rename(columns={'id_x': 'Tumor'})
rs_data = rs_data.rename(columns={'id_y': 'Normal'})


#Somatic VCF TCGA metadata
somatic_data=pd.read_csv("somaticSNP.metadata.tsv",sep='\t')
#Merge with rna-seq
somatic_data["cases.submitter_id"]=somatic_data["cases"].str[:12]
design=pd.merge(rs_data,somatic_data,on=["cases.submitter_id"])



#Generate design.txt file for RIA
design = design.rename(columns={'cases.submitter_id': 'ID'})
design=design[["Tumor","Normal","ID","id"]]
design = design.rename(columns={'id': 'VCF_ID'})
design.to_csv("design.txt",sep='\t',index=False,mode='w')
