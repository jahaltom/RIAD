###########################COSMIC VCF coding
# for var in {1..22}
# do
#    echo "$var chr$var" >> chr.txt
# done
# echo "X chrX" >> chr.txt
# echo "Y chrY" >> chr.txt
# echo "MT chrM" >> chr.txt

# bcftools index -t CosmicCodingMuts.vcf.gz
# bcftools annotate --rename-chrs chr.txt  CosmicCodingMuts.vcf.gz --output-type z -o CosmicCodingMuts.temp.vcf.gz
# mv CosmicCodingMuts.temp.vcf.gz CosmicCodingMuts.vcf.gz

# bcftools index -t -f CosmicCodingMuts.vcf.gz


import pandas as pd
import subprocess
#Read in genotypes
design=pd.read_csv("design.txt",sep='\t')
#RNA-Seq id and corresponding somatic vcf id
rs_id=design["ID"].tolist()
vcf_id=design["VCF_ID"].tolist()

path2VCF="GDCdata/TCGA-LUSC/harmonized/Simple_Nucleotide_Variation/Annotated_Somatic_Mutation/"
path2_rsVCF="output/"

for i in range(0,len(rs_id)):   
    bcftools=open(rs_id[i]+"bcftools.sh", "w")
    bcftools.write("source activate Ancestry")
    bcftools.write('\n')   
    bcftools.write("bcftools index -t "+path2_rsVCF + rs_id[i] +"_SomaticMutations/"+rs_id[i]+"_SomaticMutations.vcf.gz")
    bcftools.write('\n')
    bcftools.write("gunzip -c "+path2_rsVCF + rs_id[i] +"_SomaticMutations/"+rs_id[i]+"_SomaticMutations.vcf.gz | grep -v '#' | wc -l >> TotalRNACalls")
    bcftools.write('\n')
    bcftools.write("bcftools isec -c none  CosmicCodingMuts.vcf.gz  "+path2_rsVCF + rs_id[i] +"_SomaticMutations/"+rs_id[i]+"_SomaticMutations.vcf.gz -p dir"+rs_id[i])
    bcftools.write('\n')
    bcftools.write("grep -v '#' dir"+rs_id[i]+"/0002.vcf | wc -l >> TotalRNA.COSMIC.Calls")    
    bcftools.write('\n')
    bcftools.write("bcftools view -f 'PASS,.' " + path2_rsVCF + rs_id[i] +"_SomaticMutations/"+rs_id[i]+"_SomaticMutations.vcf.gz --output-type z -o " + path2_rsVCF + rs_id[i] +"_SomaticMutations/"+rs_id[i]+"_SomaticMutations.PASS.vcf.gz")
    bcftools.write('\n')
    bcftools.write("bcftools index -t "+path2_rsVCF + rs_id[i] +"_SomaticMutations/"+rs_id[i]+"_SomaticMutations.PASS.vcf.gz")
    bcftools.write('\n')
    bcftools.write("gunzip -c "+path2_rsVCF + rs_id[i] +"_SomaticMutations/"+rs_id[i]+"_SomaticMutations.PASS.vcf.gz | grep -v '#' | wc -l >> TotalRNA_PASS_Calls")
    bcftools.write('\n')
    bcftools.write("bcftools isec -c none  CosmicCodingMuts.vcf.gz  "+path2_rsVCF + rs_id[i] +"_SomaticMutations/"+rs_id[i]+"_SomaticMutations.PASS.vcf.gz -p dirPASS"+rs_id[i])
    bcftools.write('\n')
    bcftools.write("grep -v '#' dirPASS"+rs_id[i]+"/0002.vcf | wc -l >> TotalRNA.COSMIC._PASS_Calls")
    bcftools.write('\n')  
    bcftools.write("gunzip -c "+path2VCF + vcf_id[i] +"/*somatic_annotation.vcf.gz | grep -v '#' | wc -l >> TotalDNACalls")
    bcftools.write('\n')
    bcftools.write("bcftools isec -c none  CosmicCodingMuts.vcf.gz  "+path2VCF + vcf_id[i] +"/*somatic_annotation.vcf.gz -p dir"+vcf_id[i])
    bcftools.write('\n')
    bcftools.write("grep -v '#' dir"+vcf_id[i]+"/0002.vcf | wc -l >> TotalDNA.COSMIC.Calls")
    bcftools.write('\n')    
    bcftools.write("bcftools view -f 'PASS,.' " + path2VCF + vcf_id[i] +"/*somatic_annotation.vcf.gz --output-type z -o " + path2VCF + vcf_id[i] +"/PASS.vcf.gz")
    bcftools.write('\n')
    bcftools.write("bcftools index -t "+path2VCF + vcf_id[i] +"/PASS.vcf.gz")
    bcftools.write('\n')
    bcftools.write("gunzip -c "+path2VCF + vcf_id[i] +"/PASS.vcf.gz | grep -v '#' | wc -l >> TotalDNA_PASS_Calls")
    bcftools.write('\n')
    bcftools.write("bcftools isec -c none  CosmicCodingMuts.vcf.gz  "+path2VCF + vcf_id[i] +"/PASS.vcf.gz -p dirPASS"+vcf_id[i])
    bcftools.write('\n')
    bcftools.write("grep -v '#' dirPASS"+vcf_id[i]+"/0002.vcf | wc -l >> TotalDNA.COSMIC._PASS_Calls")
    bcftools.write('\n')   
    bcftools.write("bcftools isec -c none "+path2VCF + vcf_id[i] +"/PASS.vcf.gz "+ path2_rsVCF + rs_id[i] +"_SomaticMutations/"+rs_id[i]+"_SomaticMutations.PASS.vcf.gz  -p dirPASS"+vcf_id[i]+rs_id[i])
    bcftools.write('\n')
    bcftools.write("grep -v '#' dirPASS"+vcf_id[i]+rs_id[i]+"/0002.vcf | wc -l >> RNA_DNA_PASS_Calls")
    bcftools.write('\n')
    bcftools.write("bcftools isec -c none "+path2VCF + vcf_id[i] +"/*somatic_annotation.vcf.gz "+ path2_rsVCF + rs_id[i] +"_SomaticMutations/"+rs_id[i]+"_SomaticMutations.vcf.gz  -p dir"+vcf_id[i]+rs_id[i])
    bcftools.write('\n')
    bcftools.write("grep -v '#' dir"+vcf_id[i]+rs_id[i]+"/0002.vcf | wc -l >> RNA_DNA_Calls")
    bcftools.close()
    
    subprocess.run("bash "+rs_id[i]+"bcftools.sh", shell=True, check=True)
        
        
    
    
