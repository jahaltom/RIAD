
#Create output directory
DIR='output'

#Set path to human genome


#Read in run accession IDs from txt file. 
with open ("RAids.txt") as f:
    ra=f.read().splitlines()

	
chr_list = list(range(1, 23))
chr_list.append("X")



rule all:
	input:
		expand("{wd}/{sample}/Chr{chr}.final.vcf.gz.tbi",sample=ra, chr=chr_list,wd=DIR)
		
rule bamindex:
	input: "{wd}/{sample}/Aligned.sortedByCoord.out_star.bam"
	output: "{wd}/{sample}/Aligned.sortedByCoord.out_star.bam.bai"
	shell:
		"samtools index -b {wildcards.wd}/{wildcards.sample}/Aligned.sortedByCoord.out_star.bam"
rule var_call:
	input:
		"{wd}/{sample}/Aligned.sortedByCoord.out_star.bam.bai"
	output:
		"{wd}/{sample}/Chr{chr}.final.vcf.gz.tbi"
	shell:	
		""" 
	
		
		
		samtools view -b {wildcards.wd}/{wildcards.sample}/Aligned.sortedByCoord.out_star.bam {wildcards.chr}  > {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.bam


		gatk MarkDuplicates \
		    I= {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.bam \
		    O= {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.MarkDup.bam \
		    CREATE_INDEX= true \
		    METRICS_FILE= {wildcards.wd}/{wildcards.sample}/marked_dup_metrics.{wildcards.chr}.txt \
		    VALIDATION_STRINGENCY= SILENT
		rm {wildcards.wd}/{wildcards.sample}/*Chr{wildcards.chr}.bam

		gatk AddOrReplaceReadGroups \
		    I= {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.MarkDup.bam \
		    O= {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.Grouped.bam \
		    RGSM= sample \
		    RGLB= lib \
		    RGPL= plat \
		    RGPU= plat
		rm {wildcards.wd}/{wildcards.sample}/*Chr{wildcards.chr}.MarkDup.bam

		gatk SplitNCigarReads \
		    -I {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.Grouped.bam  \
		    -O {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.final.bam \
		    -R data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna 
		rm {wildcards.wd}/{wildcards.sample}/*Chr{wildcards.chr}.Grouped.bam



		gatk --java-options "-Xmx20g" HaplotypeCaller \
		   -I {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.final.bam \
		   -O {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.vcf.gz \
		   -R data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
		   -L data/Chr{wildcards.chr}.interval_list \
		   -ERC GVCF \
		   --native-pair-hmm-threads 4  
		rm {wildcards.wd}/{wildcards.sample}/*Chr{wildcards.chr}.final.bam


		gunzip -c {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.vcf.gz |  grep -v "0/0:0:0:0:0,0,0" > {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.final.vcf
		rm {wildcards.wd}/{wildcards.sample}/*Chr{wildcards.chr}.vcf.gz
		
		####################################
		if test -f {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.final.vcf.gz; then
		    rm {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.final.vcf.gz
		fi	
		##############################
		
	
		
		bgzip {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.final.vcf   
		bcftools index -t {wildcards.wd}/{wildcards.sample}/Chr{wildcards.chr}.final.vcf.gz
		""" 
