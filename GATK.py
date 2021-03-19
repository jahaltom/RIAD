import os
import pysam
from pyrpipe import sra,qc,mapping
from pyrpipe.runnable import Runnable


#Create output directory
working_dir='example_output'

#Set path to human genome
gen='/work/LAS/xgu-lab/Haplo/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna'
##Creating the FASTA sequence dictionary file
Dict=Runnable(command='gatk')
##specify options; these can be specified into yaml file too
param={'CreateSequenceDictionary':'', '-R': gen}
Dict.run(**param)
#Creating the fasta index file
pysam.faidx(gen)

#Create STAR index direcotry. 
star_index='star_index/HomoS'
##initialize objects
#creates a star object to use with threads.
star=mapping.Star(index=star_index,genome=gen,threads=4,**{'--twopassMode':'Basic'})
#creates a trim_galore object.
trim_galore=qc.Trimgalore()



#Pull SRR Ids from txt file. 
with open ("SRRids.txt") as f:
	SRR=f.read().splitlines()
 
#Take in 1 SRR ID at a time.     
for s in SRR:   
    #Download fastq(s) for SRR ID 
    srr_object=sra.SRA(s,directory=working_dir)
    #Run through trim_galore and STAR
    srr_object.trim(trim_galore).align(star)
    
    #Change working directory to output directory for the SRR ID. 
    os.chdir('/work/LAS/xgu-lab/Haplo/example_output/'+s)
    
    
    
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
    param={'AddOrReplaceReadGroups':'','INPUT=': s+'.DupRem.Sorted.bam','OUTPUT=': s+'.DupRem.Sorted.Grouped.bam' ,'SORT_ORDER=':'coordinate','RGSM=': s,'RGID=':'4','RGLB=':'lib1','RGPL=':'illumina','RGPU=':'unit1'}
    Group.run(**param)
    
    
    #Make a Runnable object that splits BAM file via GATK.
    Split=Runnable(command='gatk')   
    param={'SplitNCigarReads':'','-I': s+'.DupRem.Sorted.Grouped.bam','-O': s+'.Clean.bam' ,'-R': gen}
    Split.run(**param)
    
    #Index bam file before doing variant calling.
    pysam.index(s+'.Clean.bam')
    
    
    #Base Recalibration 
    # #create a Runnable object
    # RemoveDup=Runnable(command='gatk')
    # #specify orfipy options; these can be specified into orfipy.yaml too
    # param={'BaseRecalibrator':'','-I': '','-O': 'recal_data.table' ,'-R':'GCA_000001405.15_GRCh38_no_alt_analysis_set.fna'}
    # RemoveDup.run(**param)
    
    # #create a Runnable object
    # RemoveDup=Runnable(command='gatk')
    # #specify orfipy options; these can be specified into orfipy.yaml too
    # param={'SplitNCigarReads':'','-I': 'Sorted.bam','-O': 'Split.bam' ,'-R':'GCA_000001405.15_GRCh38_no_alt_analysis_set.fna'}
    # RemoveDup.run(**param)
    
    
    #Make a Runnable object that calls variants file via GATK.
    VarCall=Runnable(command='gatk')   
    param={'--java-options':'-Xmx4g', 'HaplotypeCaller':'', '-I': s+'.Clean.bam','-O': s+'.vcf.gz' ,'-R': gen}
    VarCall.run(**param)
    

    
