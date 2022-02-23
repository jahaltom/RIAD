import pandas as pd


#Create output directory
DIR='output'


#Study design with run accession IDs.
design=pd.read_csv("design",sep="\t")
normal=design['Normal'].tolist() 
tumor=design['Tumor'].tolist() 
tumor = [x for x in tumor if pd.isnull(x) == False]
ra=tumor+normal

#Get list of Uniq IDs
ids = list(set(design['ID'].tolist()))


chr_list = list(range(1, 23))



rule all:
    input: expand("{wd}/{id}.{chr}.somatic.vcf.gz",wd=DIR,chr=chr_list,id=ids)


rule BamFormating:
    input:
        expand("{wd}/{sample}/2ndPass.Aligned.sortedByCoord.out.bam",sample=ra,chr=chr_list,wd=DIR)
    output:
        "{wd}/{sample}/Chr{chr}.final.bam"

    shell:
        """
        samtools view -b {input} {wildcards.chr}  > {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.bam

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
            --known-sites data/Chr{wildcards.chr}_SNPs.vcf.gz \
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
        expand("{wd}/{sample}/Chr{chr}.final.bam",chr=chr_list,wd=DIR,normalSample=normal,sample=ra)
        
    output:
        "{wd}/{normalSample}/Chr{chr}.{normalSample}.normal.vcf.gz"
   

    shell:
        """   
        gatk Mutect2 \
          -R data/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna \
          -I {wildcards.wd}/{wildcards.normalSample}/Chr{wildcards.chr}.final.bam \
          -L {wildcards.chr} \
          --max-mnp-distance 0 \
          -O {wildcards.wd}/{wildcards.normalSample}/Chr{wildcards.chr}.{wildcards.normalSample}.normal.vcf.gz
          
          echo $"{wildcards.normalSample}\t{wildcards.wd}/{wildcards.normalSample}/Chr{wildcards.chr}.{wildcards.normalSample}.normal.vcf.gz" >> {wildcards.wd}/normals_for_pon_{wildcards.chr}_vcfs.sample_map
        """
  
rule GenomicsDBImport :
    input:
        expand("{wd}/{normalSample}/Chr{chr}.{normalSample}.normal.vcf.gz",chr=chr_list,wd=DIR,normalSample=normal)
    output:
        directory("{wd}/pon_{chr}_db")
    shell:         
        """       
         gatk GenomicsDBImport  \
           -R data/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna \
           -L {wildcards.chr} \
           --genomicsdb-workspace-path {wildcards.wd}/pon_{wildcards.chr}_db \
           --sample-name-map {wildcards.wd}/normals_for_pon_{wildcards.chr}_vcfs.sample_map
        """
        

rule CreateSomaticPanelOfNormals:
    input:
         directory(expand("{wd}/pon_{chr}_db",chr=chr_list,wd=DIR))
    output:
        "{wd}/pon.{chr}.vcf.gz"
    shell:         
        """       
         gatk CreateSomaticPanelOfNormals \
           -R data/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna \
           -V gendb://{wildcards.wd}/pon_{wildcards.chr}_db \
           -O {wildcards.wd}/pon.{wildcards.chr}.vcf.gz \
           -L {wildcards.chr}
        """

rule SomaticVarCall:
    input:      
        expand("{wd}/pon.{chr}.vcf.gz",wd=DIR,chr=chr_list,id=ids)

    output:
        "{wd}/{id}.{chr}.somatic.vcf.gz"

    run:
        
        tumorSamples=design.loc[design['ID'] == id]['Tumor'].tolist()       
        tumor_files=""        
        for i in tumorSamples:
            tumor_files= tumor_files + "-I {wildcards.wd}/" + i +"/Chr{wildcards.chr}.final.bam "
            
            
        normalSamples=design.loc[design['ID'] == id]['Normal'].tolist()
        normal_files=""
        normal_list=""        
        for i in normalSamples:
            normal_files= normal_files + "-I {wildcards.wd}/" + i +"/Chr{wildcards.chr}.final.bam "
            normal_list= normal_list + "-normal " + i + " "
            
                  
        shell("gatk Mutect2 \
          -R data/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna " + tumor_files + normal_files + normal_list +
          "-L {wildcards.chr} \
          --germline-resource data/gnomad.hg38.vcf.gz \
          -panel-of-normals {wildcards.wd}/pon.{wildcards.chr}.vcf.gz \
          --f1r2-tar-gz {wildcards.wd}/{wildcards.id}.f1r2.tar.gz \
          -O {wildcards.wd}/{wildcards.id}.{wildcards.chr}.somatic.vcf.gz")
  
          
          
          
          
          
          
          
          
          
          
  
