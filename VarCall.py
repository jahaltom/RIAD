#Create output directory
DIR='output'

#Set path to human genome




#Read in run accession IDs from txt file. 
with open ("RAids.txt") as f:
    ra=f.read().splitlines()



rule all:
	input: expand("{wd}/{sample}/vcf.gz".format(sample=s,wd=DIR) for s in ra)  
    
rule var_call:
    input:
        [directory("{wd}/{sample}").format(wd=DIR,sample=s) for s in ra]
    output: "{wd}/{sample}/vcf.gz"
    shell:
        """
        #Removes all but 1KG Chrs
        #samtools index -b {wildcards.wd}/{wildcards.sample}/Aligned.sortedByCoord.out_star.bam
        #samtools view -b {wildcards.wd}/{wildcards.sample}/Aligned.sortedByCoord.out_star.bam {{1..22}} X  > {wildcards.wd}/{wildcards.sample}/output.bam

        
        
    

        gatk MarkDuplicates \
            I= {wildcards.wd}/{wildcards.sample}/Aligned.sortedByCoord.out_star.bam \
            O= {wildcards.wd}/{wildcards.sample}/MarkDup.bam \
            CREATE_INDEX= true \
            METRICS_FILE= {wildcards.wd}/{wildcards.sample}/marked_dup_metrics.txt \
            VALIDATION_STRINGENCY= SILENT
            
        gatk AddOrReplaceReadGroups \
            I= {wildcards.wd}/{wildcards.sample}/MarkDup.bam \
            O= {wildcards.wd}/{wildcards.sample}/Grouped.bam \
            RGSM= {wildcards.sample} \
            RGLB= lib \
            RGPL= plat \
            RGPU= plat
        rm {wildcards.wd}/{wildcards.sample}/Mark*


   
        gatk SplitNCigarReads \
            -I {wildcards.wd}/{wildcards.sample}/Grouped.bam \
            -O {wildcards.wd}/{wildcards.sample}/final.bam \
            -R /work/LAS/xgu-lab/Ancestry/data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna 
          
            
        rm {wildcards.wd}/{wildcards.sample}/Grouped*
        
        gatk --java-options "-Xmx20g" HaplotypeCaller \
           -I {wildcards.wd}/{wildcards.sample}/final.bam \
           -O {wildcards.wd}/{wildcards.sample}/vcf.gz \
           -R /work/LAS/xgu-lab/Ancestry/data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
           -L /work/LAS/xgu-lab/Ancestry/data/ChrAll.interval_list \
           -ERC GVCF \
           --native-pair-hmm-threads 16  
        rm {wildcards.wd}/{wildcards.sample}/final*
        
       
        gunzip -c {wildcards.wd}/{wildcards.sample}/vcf.gz |  grep -v "0/0:0:0:0:0,0,0" > {wildcards.wd}/{wildcards.sample}/vcf
        rm {wildcards.wd}/{wildcards.sample}/vcf.gz*
        bgzip {wildcards.wd}/{wildcards.sample}/vcf        
        bcftools index -t {wildcards.wd}/{wildcards.sample}/vcf.gz
        """
#rm {wildcards.wd}/{wildcards.sample}/Aligned.sortedByCoord.out_star.bam
