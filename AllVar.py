import pandas as pd
#Create output directory
DIR='output'

#Set path to human genome


#Read in run accession IDs from txt file.
with open ("RAids.txt") as f:
    ra=f.read().splitlines()


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
    
    
    

rule all:
    input:
        expand("{wd}/{sample}/ClinVarResults.tsv",sample=ra,wd=DIR)


rule var_call:
    input:
        "{wd}/{sample}/2ndPass.Aligned.sortedByCoord.out.bam"
    output:
        "{wd}/{sample}/Chr{chr}.filtered.gvcf.gz.tbi"
    group:
        "Genotype"
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
            --known-sites data/af-only-gnomad.hg38.vcf.gz \
            -O {wildcards.wd}/{wildcards.sample}/recal_data.table
            
        gatk ApplyBQSR \
            -R data/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna \
            -I {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.split.bam \
            --bqsr-recal-file {wildcards.wd}/{wildcards.sample}/recal_data.table \
            -O {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.final.bam
        rm {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.split*
               

        gatk --java-options "-Xmx20g" HaplotypeCaller \
           -I {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.final.bam \
           -O {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.gvcf.gz \
           -R data/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna \
           -L chr{wildcards.chr} \
           -ERC GVCF \
           --native-pair-hmm-threads 7
        rm {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.final*


        gunzip -c {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.gvcf.gz |  grep -v "0/0:0:0:0:0,0,0" | bgzip > {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.filtered.gvcf.gz
        rm {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.gvcf.gz*

        bcftools index -t {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.filtered.gvcf.gz
        """
rule Genotype:
    input:
        "{wd}/{sample}/Chr{chr}.filtered.gvcf.gz.tbi"
    output:
        "{wd}/{sample}/{sample}.Chr{chr}.vcf.gz.tbi"
    group:
        "Genotype"
    shell:
        """
        gatk --java-options "-Xmx10g" GenotypeGVCFs \
        -R data/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna \
        -V {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.filtered.gvcf.gz \
        -O {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.vcf.gz \
        --include-non-variant-sites


        gatk --java-options "-Xmx10g" VariantFiltration \
        -R data/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna \
        -V {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.vcf.gz \
        -O {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.filtered.vcf.gz \
        --filter-expression "DP < 5.0" \
        --filter-name "filter_DP" \
        --filter-expression "MQ < 40.0" \
        --filter-name "filter_MQ" \
        --filter-expression "FS > 60.0" \
        --filter-name "filter_FS" \
        --filter-expression "MQRankSum < -12.5" \
        --filter-name "filter_MQRankSum" \
        --filter-expression "ReadPosRankSum < -8.0" \
        --filter-name "filter_ReadPosRankSum" \
        --filter-expression "QD < 2.0" \
        --filter-name "filter_QD"
        rm {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.vcf.gz*


        gatk --java-options "-Xmx10g" SelectVariants \
        -R data/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna \
        -V {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.filtered.vcf.gz  \
        -O {wildcards.wd}/{wildcards.sample}/{wildcards.sample}.Chr{wildcards.chr}.temp.vcf.gz \
        --exclude-filtered true
        rm {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.filtered.vcf.gz*

        gunzip -c  {wildcards.wd}/{wildcards.sample}/{wildcards.sample}.Chr{wildcards.chr}.temp.vcf.gz | grep -v  "\./\." | bgzip > {wildcards.wd}/{wildcards.sample}/{wildcards.sample}.Chr{wildcards.chr}.vcf.gz
        bcftools index -t {wildcards.wd}/{wildcards.sample}/{wildcards.sample}.Chr{wildcards.chr}.vcf.gz
        rm {wildcards.wd}/{wildcards.sample}/{wildcards.sample}.Chr{wildcards.chr}.temp.vcf.gz*
        """


rule concat:
    input:
        expand("{wd}/{sample}/{sample}.Chr{chr}.vcf.gz.tbi",sample=ra,chr=chr_list,wd=DIR)
    output:
        "{wd}/{sample}/All.{sample}.vcf.gz"
    shell:
        """
        #rm {wildcards.wd}/{wildcards.sample}/2ndPass.Aligned.sortedByCoord.out.bam*

        rm {wildcards.wd}/{wildcards.sample}/*filtered.gvcf.gz*


        bcftools concat {wildcards.wd}/{wildcards.sample}/*.vcf.gz --output-type z --threads 7 > {wildcards.wd}/{wildcards.sample}/All.{wildcards.sample}.vcf.gz
        bcftools index -t {wildcards.wd}/{wildcards.sample}/All.{wildcards.sample}.vcf.gz
        """
     
rule ClinVar:
    input: 
        "{wd}/{sample}/All.{sample}.vcf.gz"
    output: 
        "{wd}/{sample}/ClinVarResults.tsv"
    run:
        srrid=str(input).split("/")[1]
        shell("bcftools isec -c none data/clinvar.vcf.gz {wildcards.wd}/{wildcards.sample}/All.{wildcards.sample}.vcf.gz --output-type z --threads 7 -p {wildcards.wd}/{wildcards.sample}/dir")
        shell('bcftools merge {wildcards.wd}/{wildcards.sample}/dir/0002.vcf.gz {wildcards.wd}/{wildcards.sample}/dir/0003.vcf.gz | grep -v "##"   > {wildcards.wd}/{wildcards.sample}/dir/merged.vcf')
        ClinVarSummary=pd.read_csv("data/ClinVarSummaryGRCh38.txt",sep='\t')
        ClinVarSummary=ClinVarSummary.rename(columns = {'VariationID':'ID'})
        Merged=pd.read_csv(DIR+"/"+srrid + "/dir/merged.vcf",sep='\t')
        results=pd.merge(ClinVarSummary,Merged,on=['ID'])
        results.to_csv(DIR+"/"+srrid + "/ClinVarResults.tsv",mode="w", header=True,index=False,sep="\t")
        
      
     
        
        












