import pandas as pd




#Input dict
vars={
        #VariantID, Chr, Position, Ref, Alt.
        "rs1990760":	[2, 162267541 ,'C','T'],
        "rs73885319":	[22, 36265860,'A','G'],
        "rs60910145":	[22, 36265988,'T','G']
      
}


metadata=pd.read_csv("metadata.txt",sep="\t")[["run_accession","Eth1"]]
##Create dict that has self-reported ethnicity as key.
d = metadata.groupby('Eth1')['run_accession'].apply(list).to_dict()






rule all:
        input:  expand("output/{sample}/done",sample=d["POP"])
                
rule quant:
        input:  expand("output/{sample}/Chr22.{sample}.vcf.gz.tbi",sample=d["POP"])
       
        output:
            "output/{sample}/done"
        run:      
            for var in vars:
                try:
                    shell("gunzip -c output/{wildcards.sample}/Chr"+str(vars[var][0])+".{wildcards.sample}.vcf.gz |  grep -v '#' | awk '{{if($2=="+str(vars[var][1])+" && $4==\""+vars[var][2]+"\" && $5==\""+vars[var][3]+"\")print}}'  | grep  '0/1' >> Hets_POP_"   +var.replace("(", "_").replace(")", "_").replace(" ", "_").replace(">", "_"),shell=True) 
                except:
                    try: 
                        shell("gunzip -c output/{wildcards.sample}/Chr"+str(vars[var][0])+".{wildcards.sample}.vcf.gz |  grep -v '#' | awk '{{if($2=="+str(vars[var][1])+" && $4==\""+vars[var][2]+"\" && $5==\""+vars[var][3]+"\")print}}' | grep '1/1' >> HomAlt_POP_"""+var.replace("(", "_").replace(")", "_").replace(" ", "_").replace(">", "_"),shell=True)
                    except: 
                        try:
                            shell("gunzip -c output/{wildcards.sample}/Chr"+str(vars[var][0])+".{wildcards.sample}.vcf.gz |  grep -v '#' | awk '{{if($2=="+str(vars[var][1])+" && $4==\""+vars[var][2]+"\" && $5==\""+vars[var][3]+"\")print}}' | grep  '0/0' >> HomRef_POP_""" +var.replace("(", "_").replace(")", "_").replace(" ", "_").replace(">", "_"),shell=True)
                        except:
                            pass 
