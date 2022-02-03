import pandas as pd


#Create output directory
DIR='output'
somatic_calls='Somatic_Calls'

#Study design with run accession IDs.
design=pd.read_csv("design",sep="\t")
normal=design['Normal'].tolist() 
tumor=design['Tumor'].tolist() 
ra=tumor+normal

chr_list = list(range(1, 23))
chr_list.append("M")




rule BamFormating:
    input:
        expand("{wd}/{sample}/2ndPass.Aligned.sortedByCoord.out.bam",sample=ra,chr=chr_list,wd=DIR)
    output:
        "{wd}/{sample}/Chr{chr}.filtered.gvcf.gz.tbi"
    group:
        "SomaticVarCall"
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
            -O {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.final.bam \
            -R data/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna
        rm {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.Grouped.bam
        
        gatk BaseRecalibrator \
            -I {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.split.bam \
            -R data/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna \
            --known-sites data/Chr{wildcards.chr}_SNPs.vcf.gz \
            -O {wildcards.wd}/{wildcards.sample}/recal_data.table
            
        gatk ApplyBQSR \
            -R data/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna \
            -I {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.split.bam \
            --bqsr-recal-file {wildcards.wd}/{wildcards.sample}/recal_data.table \
            -O {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.final.bam
        rm {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.split.bam
        """


rule CreateSomaticPanelOfNormals:
    input:
        expand("{wd}/{normalSample}/Chr{wildcards.chr}.final.bam",normalSample=normal,chr=chr_list,wd=DIR,somatic_calls=somatic_calls)
    output:
        
    group:
        "SomaticVarCall"
    shell:
        """
        gatk Mutect2 \
          -R data/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna \
          -I {wildcards.wd}/{wildcards.normalSample}/Chr1.final.bam \
          -tumor {wildcards.normalSample} \
          -L {wildcards.chr} \
          --germline-resource data/Chr{wildcards.chr}_SNPs.vcf.gz \
          -O {wildcards.wd}/{wildcards.normalSample}/{wildcards.normalSample}.normal.vcf.gz
          
          
         echo "{wildcards.wd}/{wildcards.normalSample}/{wildcards.normalSample}.normal.vcf.gz" >> {wildcards.wd}/{wildcards.somatic_calls}/normals_for_pon_vcf.args
         
         gatk CreateSomaticPanelOfNormals \
           -vcfs {wildcards.wd}/{wildcards.somatic_calls}/normals_for_pon_vcf.args \
           -O {wildcards.wd}/{wildcards.somatic_calls}/pon.vcf.gz
        """

rule SomaticVarCall:
    input:
        expand("{wd}/{tumorSample}/Chr{wildcards.chr}.final.bam",tumorSample=tumor,chr=chr_list,wd=DIR,somatic_calls=somatic_calls),
    output:
        
    group:
        "SomaticVarCall"
    shell:
        """
        gatk Mutect2 \
          -R data/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna \
          -I {wildcards.wd}/{wildcards.tumorSample}/Chr1.final.bam \
          -L {wildcards.chr}
          --germline-resource data/Chr{wildcards.chr}_SNPs.vcf.gz \
          -panel-of-normals {wildcards.wd}/{wildcards.somatic_calls}/pon.vcf.gz \
          -O {wildcards.wd}/{wildcards.somatic_calls}/{wildcards.tumorSample}.somatic.vcf.gz
        """
