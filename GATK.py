import os
import pysam
from pyrpipe import sra,qc,mapping
from pyrpipe.runnable import Runnable
import pandas as pd

#snakemake -j 240 --cluster "sbatch -t 05:00:00 -c 20 -p RM-shared"


#Create output directory
DIR='output'

#Set path to human genome

gen='/ocean/projects/mcb200036p/jahaltom/Ancestry/data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna'


interval='/ocean/projects/mcb200036p/jahaltom/Ancestry/data/ChrAllIntervals.interval_list'

chr_list = list(range(1, 23))
chr_list.append("X")


#Read in run accession IDs from txt file. 
with open ("RAids.txt") as f:
    ra=f.read().splitlines()

wildcard_constraints:
    chr= '|'.join([str(x) for x in chr_list])
            	

rule all:
	input: expand("{wd}/{sample}/vcf.gz".format(sample=s,wd=DIR) for s in ra)
  
    
rule var_call:
    input:
        ["{wd}/{sample}/Aligned.sortedByCoord.out_star.bam".format(wd=DIR,sample=s) for s in ra]
    output: "{wd}/{sample}/vcf.gz"
    run:
        srrid=str(output).split("/")[1]
        path=DIR + "/" + srrid + "/"
      
     
        #Make a Runnable object that removes duplicates via picard. 
        RemoveDup=Runnable(command='gatk')   
        param={'MarkDuplicates':'','INPUT=': path + 'Aligned.sortedByCoord.out_star.bam','OUTPUT=': path + 'DupRem.bam' ,'CREATE_INDEX=':'true','METRICS_FILE=': path + 'marked_dup_metrics.txt','VALIDATION_STRINGENCY=':'SILENT'}
        RemoveDup.run(**param)
        #shell("rm {wildcards.wd}/{wildcards.sample}/Aligned.sortedByCoord.out_star.bam")
        
        
        #Make a Runnable object that adds read group info to BAM file and sorts it via picard.
        Group=Runnable(command='gatk')  
        param={'AddOrReplaceReadGroups':'','INPUT=': path + 'DupRem.bam','OUTPUT=': path + 'DupRem.Sorted.Grouped.bam' ,'SORT_ORDER=':'coordinate','RGSM=': srrid,'RGID=':'4','RGLB=':'lib1','RGPL=':'illumina','RGPU=':'unit1'}
        Group.run(**param)
        shell("rm {wildcards.wd}/{wildcards.sample}/DupRem.bam")
        
        #Make a Runnable object that splits BAM file via GATK.
        Split=Runnable(command='gatk')   
        param={'SplitNCigarReads':'','-I': path + 'DupRem.Sorted.Grouped.bam','-O': path + 'final.bam' ,'-R': gen}
        Split.run(**param)
        shell("rm {wildcards.wd}/{wildcards.sample}/DupRem.Sorted.Grouped.bam")
        #Index bam file before doing variant calling.
        pysam.index(path+ 'final.bam')
                
        
                #Make a Runnable object that calls variants file via GATK.
        VarCall=Runnable(command='gatk')
        param={'--java-options':'-Xmx20g', 'HaplotypeCaller':'', '-I': path + 'final.bam','-O': path + 'vcf.gz' ,'-R': gen,'-L': interval,'-ERC':'GVCF','--native-pair-hmm-threads': '10'}
        VarCall.run(**param)
        shell("rm {wildcards.wd}/{wildcards.sample}/final.bam")
        # shell("gunzip -c vcf.gz |  grep -v ""0/0:0:0:0:0,0,0"" > vcf2.gz")
        # shell("mv vcf2.gz vcf.gz")
        # shell("rm vcf.gz.tbi")
        # shell("bcftools index -t vcf.gz ")

# rule JointGenotype:
#     input: ["Ancestry/{wd}/{sample}/vcf.gz".format(wd=DIR,sample=s) for s in ra]
#     output: "joint/SNP_filtered.chr{chr}.joint.vcf.gz"
    
#     shell:
#         """  
#         gatk --java-options "-Xmx200g -Xms200g" \
#          GenomicsDBImport \
#         --genomicsdb-workspace-path joint/chr{wildcards.chr}_database \
#         --batch-size 2 \
#      	-L {wildcards.chr} \
#         --sample-name-map Ancestry/data/cohort.sample_map \
#         --tmp-dir temp \
#         --reader-threads 40
        
#         gatk --java-options "-Xmx200g" GenotypeGVCFs \
#         -R /ocean/projects/mcb200036p/jahaltom/Ancestry/data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
#         -V gendb://joint/chr{wildcards.chr}_database \
#         -O joint/chr{wildcards.chr}.joint.vcf.gz \
#         -L {wildcards.chr} \
#         --include-non-variant-sites \
#         --tmp-dir temp 
#         rm -r joint/chr{wildcards.chr}_database
                
#         gatk --java-options "-Xmx100g" VariantFiltration \
#         -R /ocean/projects/mcb200036p/jahaltom/Ancestry/data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
#         -V joint/chr{wildcards.chr}.joint.vcf.gz \
#         -O joint/filtered.chr{wildcards.chr}.joint.vcf.gz \
#         --filter-expression "DP < 10.0" \
#         --filter-name "filter_DP" \
#         --filter-expression "MQ < 40.0" \
#         --filter-name "filter_MQ" \
#         --filter-expression "FS > 60.0" \
#         --filter-name "filter_FS" \
#         --filter-expression "MQRankSum < -12.5" \
#         --filter-name "filter_MQRankSum" \
#         --filter-expression "ReadPosRankSum < -8.0" \
#         --filter-name "filter_ReadPosRankSum" \
#         --filter-expression "QD < 2.0" \
#         --filter-name "filter_QD"
#         rm joint/chr{wildcards.chr}.joint.vcf.gz*
        
        
#         gatk --java-options "-Xmx100g" SelectVariants \
#         -R /ocean/projects/mcb200036p/jahaltom/Ancestry/data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
#         -V joint/filtered.chr{wildcards.chr}.joint.vcf.gz \
#         -O joint/SNP_filtered.chr{wildcards.chr}.joint.vcf.gz \
#         --exclude-filtered true \
#         --select-type-to-exclude INDEL      
        
#         rm joint/filtered.chr{wildcards.chr}.joint.vcf.gz*
     
#         for file in joint/SNP_filtered.chr{wildcards.chr}.joint.vcf.gz; do
#           for sample in `bcftools query -l $file`; do
#             bcftools view -Oz -s $sample -o $sample.{wildcards.chr}.vcf.gz $file
#           done
#         done
#         rm joint/SNP_filtered.chr{wildcards.chr}.joint.vcf.gz*
                    
#         """
        

            


	
