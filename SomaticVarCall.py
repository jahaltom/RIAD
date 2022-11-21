import pandas as pd
import os



#Read in config file
configfile: "config.yaml"
OutputDir=config['OutputDir']
Interval = config['Interval']
bcftools_threads=config['bcftools_threads']

chr_list=[]
if Interval == "All":
    for i in range(1, 23):
        chr_list.append(str(i))
elif Interval == "AllMYX":
    for i in range(1, 23):
        chr_list.append(str(i))
    chr_list.append("M")
    chr_list.append("Y")
    chr_list.append("X")     
else:
    chr_list=Interval
    
    


#Study design with run accession IDs.
design=pd.read_csv("design.txt",sep="\t")
normal=design['Normal'].tolist() 
tumor=design['Tumor'].tolist() 
tumor = [x for x in tumor if pd.isnull(x) == False]
ra=tumor+normal

#Get list of Uniq IDs
ids = list(set(design['ID'].tolist()))









rule all:
    input: expand("{wd}/{id}_SomaticMutations/{id}_SomaticMutations.vcf.gz",wd=OutputDir,id=ids)


rule BamFormating:
    input:
        expand("{wd}/{sample}/2ndPass.Aligned.sortedByCoord.out.bam",sample=ra,chr=chr_list,wd=OutputDir)
    output:
        "{wd}/{sample}/Chr{chr}.final.bam"

    shell:
        """
        source activate Ancestry
        samtools view -b {wildcards.wd}/{wildcards.sample}/2ndPass.Aligned.sortedByCoord.out.bam chr{wildcards.chr}  > {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.bam

        gatk MarkDuplicates \
            I= {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.bam \
            O= {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.MarkDup.bam \
            CREATE_INDEX= true \
            METRICS_FILE= {wildcards.wd}/{wildcards.sample}/marked_dup_metrics.{wildcards.chr}.txt \
            VALIDATION_STRINGENCY= SILENT
        rm {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.bam

        gatk AddOrReplaceReadGroups \
            I= {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.MarkDup.bam \
            O= {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.Grouped.bam \
            RGSM= {wildcards.sample} \
            RGLB= lib \
            RGPL= plat \
            RGPU= plat
        rm {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.MarkDup*

        gatk SplitNCigarReads \
            -I {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.Grouped.bam  \
            -O {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.split.bam \
            -R data/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna
        rm {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.Grouped.bam
        
        gatk BaseRecalibrator \
            -I {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.split.bam \
            -R data/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna \
            --known-sites data/af-only-gnomad.hg38.vcf.gz \
            -O {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.recal_data.table
            
        gatk ApplyBQSR \
            -R data/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna \
            -I {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.split.bam \
            --bqsr-recal-file {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.recal_data.table \
            -O {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.final.bam
        rm {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.split*
        """


rule VarCallNormals:
    input:
        expand("{wd}/{sample}/Chr{chr}.final.bam",chr=chr_list,wd=OutputDir,normalSample=normal,sample=ra)
        
    output:
        "{wd}/{normalSample}/Chr{chr}.{normalSample}.normal.vcf.gz"
   

    shell:
        """   
        source activate Ancestry
        gatk Mutect2 \
          -R data/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna \
          -I {wildcards.wd}/{wildcards.normalSample}/Chr{wildcards.chr}.final.bam \
          -L chr{wildcards.chr} \
          --max-mnp-distance 0 \
          -O {wildcards.wd}/{wildcards.normalSample}/Chr{wildcards.chr}.{wildcards.normalSample}.normal.vcf.gz
          
          echo $"{wildcards.normalSample}\t{wildcards.wd}/{wildcards.normalSample}/Chr{wildcards.chr}.{wildcards.normalSample}.normal.vcf.gz" >> {wildcards.wd}/normals_for_pon_{wildcards.chr}_vcfs.sample_map
        """
  
rule GenomicsDBImport :
    input:
        expand("{wd}/{normalSample}/Chr{chr}.{normalSample}.normal.vcf.gz",chr=chr_list,wd=OutputDir,normalSample=normal)
    output:
        directory("{wd}/pon_{chr}_db")
    shell:         
        """       
         source activate Ancestry
         gatk GenomicsDBImport  \
           -R data/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna \
           -L chr{wildcards.chr} \
           --genomicsdb-workspace-path {wildcards.wd}/pon_{wildcards.chr}_db \
           --sample-name-map {wildcards.wd}/normals_for_pon_{wildcards.chr}_vcfs.sample_map
        """
        

rule CreateSomaticPanelOfNormals:
    input:
         directory(expand("{wd}/pon_{chr}_db",chr=chr_list,wd=OutputDir))
    output:
        "{wd}/pon.{chr}.vcf.gz"
    shell:         
        """
         source activate Ancestry
         gatk CreateSomaticPanelOfNormals \
           -R data/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna \
           -V gendb://{wildcards.wd}/pon_{wildcards.chr}_db \
           -O {wildcards.wd}/pon.{wildcards.chr}.vcf.gz \
          --germline-resource data/af-only-gnomad.hg38.vcf.gz \
           -L chr{wildcards.chr}
        """


##Joint calling tumor and normal samples from the same individual.
rule SomaticVarCall:
    input:
        expand("{wd}/pon.{chr}.vcf.gz",wd=OutputDir,chr=chr_list,id=ids)

    output:
        "{wd}/{id}_SomaticMutations/{id}.Chr.{chr}.somatic.unfiltered.vcf.gz"

    run:

        try:
            os.mkdir(OutputDir+"/"+ wildcards.id + "_SomaticMutations")
        except: pass

        tumorSamples=design.loc[design['ID'] == wildcards.id]['Tumor'].tolist()
        tumor_files=""
        for i in tumorSamples:
            tumor_files= tumor_files + "-I "+wildcards.wd+"/" + i +"/Chr"+wildcards.chr+".final.bam "


        normalSamples=design.loc[design['ID'] == wildcards.id]['Normal'].tolist()
        normal_files=""
        normal_list=""
        for i in normalSamples:
            normal_files= normal_files + "-I "+wildcards.wd+"/" + i +"/Chr"+wildcards.chr+".final.bam "
            normal_list= normal_list + "-normal " + i + " "








        gatk=open("gatk."+wildcards.chr+"."+wildcards.id+".sh", "w")
        gatk.write("source activate Ancestry")
        gatk.write('\n')
        gatk.write("gatk Mutect2 \
          -R data/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna " + tumor_files + normal_files + normal_list +
          "-L chr"+wildcards.chr+" \
          --germline-resource data/af-only-gnomad.hg38.vcf.gz \
          -panel-of-normals "+wildcards.wd+"/pon."+wildcards.chr+".vcf.gz \
          --f1r2-tar-gz "+wildcards.wd+"/"+wildcards.id+"_SomaticMutations/"+wildcards.id+".Chr."+wildcards.chr+".f1r2.tar.gz \
          -O "+wildcards.wd+"/"+wildcards.id+"_SomaticMutations/"+wildcards.id+".Chr."+wildcards.chr+".somatic.unfiltered.vcf.gz")
        gatk.close()
        subprocess.run("bash gatk."+wildcards.chr+"."+wildcards.id+".sh", shell=True, check=True)



rule CalculateContamination:
    input:
        expand("{wd}/{id}_SomaticMutations/{id}.Chr.{chr}.somatic.unfiltered.vcf.gz",wd=OutputDir,chr=chr_list,tumorSamples=tumor,id=ids)
    output:
        "{wd}/{tumorSamples}/Chr{chr}.calculatecontamination.table"

    run:

        idx=design.index[design['Tumor']==wildcards.tumorSamples].tolist()
        normal_sample=design["Normal"][idx][idx[0]]


        gatk=open("gatk."+wildcards.chr+"."+wildcards.tumorSamples+".sh", "w")
        gatk.write("source activate Ancestry")
        gatk.write('\n')
        gatk.write("gatk GetPileupSummaries \
        -I "+wildcards.wd+"/"+wildcards.tumorSamples+"/Chr"+wildcards.chr+".final.bam \
        -V data/af-only-gnomad.hg38.vcf.gz \
        -L chr"+wildcards.chr+" \
        -O "+wildcards.wd+"/"+wildcards.tumorSamples+"/getpileupsummaries.Chr"+wildcards.chr+".table")
        gatk.write('\n')
        gatk.write("gatk GetPileupSummaries \
        -I "+wildcards.wd+"/"+normal_sample+"/Chr"+wildcards.chr+".final.bam \
        -V data/af-only-gnomad.hg38.vcf.gz \
        -L chr"+wildcards.chr+" \
        -O "+wildcards.wd+"/"+normal_sample+"/getpileupsummaries.Chr"+wildcards.chr+".table")
        gatk.write('\n')
        gatk.write("gatk CalculateContamination \
        -I "+wildcards.wd+"/"+wildcards.tumorSamples+"/getpileupsummaries.Chr"+wildcards.chr+".table \
        -matched "+wildcards.wd+"/"+normal_sample+"/getpileupsummaries.Chr"+wildcards.chr+".table \
        -tumor-segmentation "+wildcards.wd+"/"+wildcards.tumorSamples+"/Chr"+wildcards.chr+".segments.table \
        -O "+wildcards.wd+"/"+wildcards.tumorSamples+"/Chr"+wildcards.chr+".calculatecontamination.table")
        gatk.close()
        subprocess.run("bash gatk."+wildcards.chr+"."+wildcards.tumorSamples+".sh", shell=True, check=True)




rule FilterCalls:
    input:
        expand("{wd}/{tumorSample}/Chr{chr}.calculatecontamination.table",wd=OutputDir,chr=chr_list,id=ids,tumorSample=tumor)
    output:
        "{wd}/{id}_SomaticMutations/{id}.Chr.{chr}.somatic.filtered.vcf.gz"
    run:

        tumorSamples=design.loc[design['ID'] == wildcards.id]['Tumor'].tolist()
        tumor_Segfiles=""
        tumor_Confiles=""
        for i in tumorSamples:
            tumor_Segfiles= tumor_Segfiles + "--tumor-segmentation "+wildcards.wd+"/" + i +"/Chr"+wildcards.chr+".segments.table "
            tumor_Confiles= tumor_Confiles + "--contamination-table "+wildcards.wd+"/" + i +"/Chr"+wildcards.chr+".calculatecontamination.table "



        gatk=open("gatk."+wildcards.chr+"."+wildcards.id+".sh", "w")
        gatk.write("source activate Ancestry")
        gatk.write('\n')
        gatk.write("gatk LearnReadOrientationModel \
          -I "+wildcards.wd+"/"+wildcards.id+"_SomaticMutations/"+wildcards.id+".Chr."+wildcards.chr+".f1r2.tar.gz \
          -O "+wildcards.wd+"/"+wildcards.id+"_SomaticMutations/read-orientation-model.Chr"+wildcards.chr+".tar.gz")
        gatk.write('\n')
        gatk.write("gatk FilterMutectCalls \
          -R data/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna \
          -V "+wildcards.wd+"/"+wildcards.id+"_SomaticMutations/"+wildcards.id+".Chr."+wildcards.chr+".somatic.unfiltered.vcf.gz \
          --stats "+wildcards.wd+"/"+wildcards.id+"_SomaticMutations/"+wildcards.id+".Chr."+wildcards.chr+".somatic.unfiltered.vcf.gz.stats \
          -L chr"+wildcards.chr+" \
          --ob-priors "+wildcards.wd+"/"+wildcards.id+"_SomaticMutations/read-orientation-model.Chr"+wildcards.chr+".tar.gz " + tumor_Segfiles + tumor_Confiles +
          "-O "+wildcards.wd+"/"+wildcards.id+"_SomaticMutations/"+wildcards.id+".Chr."+wildcards.chr+".somatic.filtered.vcf.gz")
        gatk.close()
        subprocess.run("bash gatk."+wildcards.chr+"."+wildcards.id+".sh", shell=True, check=True)

               
          
rule concat:
    input: 
        expand("{wd}/{id}_SomaticMutations/{id}.Chr.{chr}.somatic.filtered.vcf.gz",wd=OutputDir,id=ids,chr=chr_list)
    output: 
        "{wd}/{id}_SomaticMutations/{id}_SomaticMutations.vcf.gz"
    shell:
        """
        source activate Ancestry
        bcftools concat {wildcards.wd}/{wildcards.id}_SomaticMutations/*somatic.filtered.vcf.gz --output-type z --threads {config[bcftools_threads]} > {wildcards.wd}/{wildcards.id}_SomaticMutations/{wildcards.id}_SomaticMutations.vcf.gz
        bcftools index -t {wildcards.wd}/{wildcards.id}_SomaticMutations/{wildcards.id}_SomaticMutations.vcf.gz
        """
         
          
          
          
          
  
