from pyrpipe import sra,qc
import glob



#Read in config file
configfile: "config.yaml"
DIR = config['OutputDir']
THREADS=config['THREADS']
FileType=config['FileType']




#creates a trim_galore object.
trim_galore=qc.Trimgalore(threads=4)

def STAR():
    #Delet origonal fastq(s)
    shell("ls {wildcards.wd}/{wildcards.sample}/*fastq | cat | grep -v trim | xargs rm")
    #STAR 1st pass alignment
    shell(" STAR --runThreadN "+THREADS+" \
    --genomeDir data/star_index_Human \
    --outSAMtype None \
    --limitSjdbInsertNsj 5041695 \
    --outFileNamePrefix {wildcards.wd}/{wildcards.sample}/1stPass. \
    --readFilesIn {wildcards.wd}/{wildcards.sample}/*trimgalore.fastq")
    
    
#Read in run_accession or smaple ids from txt file.
with open ("ids.txt") as f:
    id=f.read().splitlines()


rule all:
        input: expand("{wd}/{sample}/2ndPass.Aligned.sortedByCoord.out.bam",sample=id,wd=DIR)

    
rule STAR1st_pass:
    output:
        "{wd}/{sample}/1stPass.SJ.out.tab"
    run:
      
        #Pathway for fastq files and salmon output 
        path=wildcards.wd + "/" + wildcards.sample+ "/"
                    
        if (FileType == 'fastq' or FileType == 'fq'):
            
            #fastq file paths
            my_files_path = glob.glob(path+'*.'+FileType)

            if len(my_files_path) == 1:
                #Run Salmon on sra object(fastq files) and delete fastq when finished.
                sra.SRA(fastq=my_files_path[0],directory=path).trim(trim_galore)
                STAR()

            elif len(my_files_path) == 2:
                #Run Salmon on sra object(fastq files) and delete fastq when finished.
                sra.SRA(fastq=my_files_path[0],fastq2=my_files_path[1],directory=path).trim(trim_galore)
                STAR()     
                        
        elif (FileType == 'SRA'):
            #Download fastq(s) from SRA. Run through trim_galore
            sra.SRA(wildcards.sample,directory=wildcards.wd).trim(trim_galore)
            STAR()
  
   
rule SJs:
    input: expand("{wd}/{sample}/1stPass.SJ.out.tab",sample=id,wd=DIR)
    output: "{wd}/all.SJ.out.tab"
    shell:
         #Combines all SJ.out.tab files from the first pass and removes SJ's that are supported by 2 or less unique mappers. Also removes annotated jusntions since these are already recorded in the genome build.
         "awk -f  data/sjCollapseSamples.awk  {input} | sort -k1,1V -k2,2n -k3,3n | awk '{{if($7>2 && $6==0)print}}' > {wildcards.wd}/all.SJ.out.tab"


rule STAR2nd_pass:
    input: "{wd}/all.SJ.out.tab"
    output:
        "{wd}/{sample}/2ndPass.Aligned.sortedByCoord.out.bam"
    run:
            #STAR 2nd pass alignment
            shell(" STAR --runThreadN "+THREADS+" \
            --genomeDir data/star_index_Human \
            --sjdbFileChrStartEnd {wildcards.wd}/all.SJ.out.tab \
            --limitSjdbInsertNsj 5041695 \
            --outSAMtype BAM SortedByCoordinate \
            --readFilesIn {wildcards.wd}/{wildcards.sample}/*trimgalore.fastq \
            --outFileNamePrefix {wildcards.wd}/{wildcards.sample}/2ndPass.")
            shell("rm {wildcards.wd}/{wildcards.sample}/*fastq")
            shell("samtools index -b {output}")
   
