import pandas as pd


# cat POP | while read i;do
# 	sed "s/POP/$i/g" snakefile  > snakefile$i 
# 	snakemake -j 300 -k -s snakefile$i --cluster "sbatch -t 01:00:00 -c 10 -p RM-shared"
# done



#https://clinvarminer.genetics.utah.edu/variants-by-gene/APOL1/significance/risk%20factor
#https://www.ebi.ac.uk/gwas/variants/rs6102095  https://www.nature.com/articles/s41431-019-0483-5
#NM_000518.5(HBB):c.20A>T (p.Glu7Val) https://www.ncbi.nlm.nih.gov/clinvar/variation/15333/?new_evidence=false
# rs1990760  https://www.mdpi.com/2227-9059/10/3/549/htm https://elifesciences.org/articles/73012  https://www.ncbi.nlm.nih.gov/snp/rs1990760

#APOL1
vars={
      #  "NM_000518.5(HBB):c.20A>T (p.Glu7Val)": [11, 5227002 , 'A', 'T'],
        "rs1990760":	[2, 162267541 ,'C','T'],
     #   "rs6102095": [20, 40692111, 'G','A'],
        "rs73885319":	[22, 36265860,'A','G'],
        "rs60910145":	[22, 36265988,'T','G']
      
}


metadata=pd.read_csv("metadata",sep="\t")[["run_accession","Eth1"]]
##Create dict that has Condition as key.
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
