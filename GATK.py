import os
import pysam
from pyrpipe import sra,qc,mapping
from pyrpipe.runnable import Runnable

#SRR11830082

#Create output directory
output_dir='output'

#Set path to human genome
gen='/work/LAS/xgu-lab/Haplo/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna'
##Creating the FASTA sequence dictionary file
#Dict=Runnable(command='gatk')
##specify options; these can be specified into yaml file too
#param={'CreateSequenceDictionary':'', '-R': gen}
#Dict.run(**param)
#Creating the fasta index file
#pysam.faidx(gen)

#Create STAR index direcotry. 
star_index='star_index/HomoS'
##initialize objects
#creates a star object to use with threads.
star=mapping.Star(index=star_index,genome=gen,threads=4,**{'--twopassMode':'Basic'})
#creates a trim_galore object.
trim_galore=qc.Trimgalore()



#Read in run accession IDs from txt file. 
with open ("SRRids.txt") as f:
	SRR=f.read().splitlines()
 
#Take in 1 run accession ID at a time.     
for s in SRR:   
	#Download fastq(s) for run accession ID 
	srr_object=sra.SRA(s,directory=output_dir)
	#Run through trim_galore and STAR
	srr_object.trim(trim_galore).align(star)

	#Change working directory to output directory for the SRR ID. 
	os.chdir('/work/LAS/xgu-lab/Haplo/output/'+s)



	#Make a Runnable object that removes duplicates via picard. 
	RemoveDup=Runnable(command='picard')   
	param={'MarkDuplicates':'','INPUT=': 'Aligned.sortedByCoord.out_star.bam','OUTPUT=': s+'.DupRem.bam' ,'REMOVE_DUPLICATES=':'true','METRICS_FILE=': 'marked_dup_metrics.txt'}
	RemoveDup.run(**param)



	# #Make a Runnable object that sorts BAM file via picard.
	# Sort=Runnable(command='picard')  
	# param={'SortSam':'','INPUT=': s+'.DupRem.bam','OUTPUT=': s+'.DupRem.Sorted.bam' ,'SORT_ORDER=':'coordinate'}
	# Sort.run(**param)

	#Make a Runnable object that adds read group info to BAM file and sorts it via picard.
	Group=Runnable(command='picard')  
	param={'AddOrReplaceReadGroups':'','INPUT=': s+'.DupRem.bam','OUTPUT=': s+'.DupRem.Sorted.Grouped.bam' ,'SORT_ORDER=':'coordinate','RGSM=': s,'RGID=':'4','RGLB=':'lib1','RGPL=':'illumina','RGPU=':'unit1'}
	Group.run(**param)


	#Make a Runnable object that splits BAM file via GATK.
	Split=Runnable(command='gatk')   
	param={'SplitNCigarReads':'','-I': s+'.DupRem.Sorted.Grouped.bam','-O': s+'.Clean.bam' ,'-R': gen}
	Split.run(**param)

	#Index bam file before doing variant calling.
	pysam.index(s+'.Clean.bam')


	#Base Recalibration 
	##create a Runnable object
	#RemoveDup=Runnable(command='gatk')    
	#param={'BaseRecalibrator':'','-I': '','-O': 'recal_data.table' ,'-R':'GCA_000001405.15_GRCh38_no_alt_analysis_set.fna'}
	# RemoveDup.run(**param)




	#Make a Runnable object that calls variants file via GATK.
	VarCall=Runnable(command='gatk')
	param={'--java-options':'-Xmx4g', 'HaplotypeCaller':'', '-I': s+ '.Clean.bam','-O': 'SRR12850399.vcf.gz' ,'-R': gen,'-L':'chr1.vcf','-ERC':'GVCF'}
	VarCall.run(**param)

	GVCF=Runnable(command='gatk')
	param={'--java-options':'-Xmx4g', 'GenotypeGVCFs':'', '-V': 'SRR12850399.vcf.gz','-O': 'SRR12850399.gvcf.vcf.gz' ,'-R': gen,'--include-non-variant-sites': '' }
	GVCF.run(**param)


VariantFiltration =Runnable(command='gatk')
param={'--java-options':'-Xmx4g', 'VariantFiltration':'', '-V': 'SRR12850399.gvcf.vcf.gz','-O': 'SRR12850399.finshed.vcf.gz' ,'-R': gen,'--filter-name':'rawwww' ,'--filter-expression': 'DP' + '<' + '10'  }
VariantFiltration .run(**param)


SelectVariants  =Runnable(command='gatk')
param={'--java-options':'-Xmx4g', 'SelectVariants':'', '-V': 'SRR12850399.finshed.vcf.gz','-O': 'SRR12850399.DONE4.vcf.gz' ,'-R': gen,'--exclude-filtered': 'true', '--select-type-to-exclude': 'INDEL' }
SelectVariants  .run(**param)
#gunzip -c SRR12850399.DONE4.vcf.gz | grep -v  '\./\.' > SRR12850399.filtered.vcf
#bgzip SRR12850399.filtered.vcf
#tabix -fp vcf SRR12850399.filtered.vcf.gz

#Must remove indels from main chr rawww
## bcftools isec -c all chr1.filtered.vcf.gz.recode.vcf.gz SRR12850399.filtered.vcf.gz -p dir
	
	
    #VariantFiltration =Runnable(command='gatk')
    #param={'--java-options':'-Xmx4g', 'VariantFiltration':'', '-V': 'SRR12850399.done.vcf.gz','-O': 'SRR12850399.finshed.vcf.gz' ,'-R': gen,'--filter-name':'rawwww' ,'--filter-expression': 'ReadPosRankSum > -8 || QD > 5 || DP > 10 ||  FS < 60 || MQ > 40 || MQRankSum > -12.5' }
    #VariantFiltration .run(**param)

	
