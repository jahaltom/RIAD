
# RNA-Seq Inferred Ancestry (RIA)
RIA is a method for infering super-population (Africa, Europe, South Asia, East Asia, and America) identity from Human RNA-seq data.
RIA leverages data from 1000 genomes project and utilizes a machine learning approach that involves principal component analysis and support vector machine. 


![alt text](https://github.com/jahaltom/RNA-Seq-Ancestry-Inference/blob/main/FlowChart.png?raw=true)

## Prerequisites

Using Conda 4.10.3, create the conda enviroment and activate:
```
conda env create -f environment.yml
conda activate Ancestry
```
or you can use the Singularity image:


```
singularity pull ria.sif library://aseetharam/ancestry/ria:latest
```

you can access the tools inside the container by prefixing:

```
module load singularity
singularity exec --bind $PWD ria.sif snakemake 
```

## Data Preparation

**1000 Genomes Project:**
The snakemake script "Prepare_1KGP" downloads chr(1-22) level VCF files from 1000 Genomes Project phase 3 on GRCh38 (https://www.internationalgenome.org/data-portal/data-collection/grch38, https://doi.org/10.12688/wellcomeopenres.15126.2) while filtering out indels. It also indexes and creates a BED for each filtered VCF file. 
```
snakemake -j 22 -s Prepare_1KGP --cluster "sbatch -t 01:00:00 -c 4 -N 1"
```

**GRCh38 Reference Genome**
The bash script "Prepare_Reference_Genome" will download the Human genome GRCh38 fasta(GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz) and the corresponding gtf, and will create a seqence dictionary and index file for the fasta. It also creates a STAR index.
```
sbatch Prepare_Reference_Genome
```

## STAR.py: 
This snakemake script performs STAR 2-pass alignment when given fastq file(s) (fastq mode) or run accession ID(s) (SRA mode).  In SRA mode, this script fetches the raw fastq files from the SRA and then uses Trimgalore for QC. Trimgalore is used for QC in fastq mode as well. After QC, the reads are then ran through a STAR 2-Pass alignment for enhanced novel SJ detection. The SJ.out.tab file for the 2nd pass is made by combining all SJ.out.tab files from the first pass and removing SJ's that are supported by 2 or less unique mappers. 

* User must specify FileType, OutputDir, and STAR_THREADS in config.yaml.
* For SRA mode, ids.txt must contain run accession ID(s).

ids.txt
```
SRR975601
SRR975602
SRR975603
SRR975604
```

* For fastq mode ids.txt must contain directory names for individual sample fastq(s). Below would be (Sample1,Sample2)

#Single-end
```
OutputDir/Sample1/Sample1.fastq
OutputDir/Sample2/Sample2.fastq
```
#Paired-end
```
OutputDir/Sample1/Sample1.r1.fastq
OutputDir/Sample1/Sample1.r2.fastq
OutputDir/Sample2/Sample2.r1.fastq
OutputDir/Sample2/Sample2.r2.fastq
```

ids.txt
```
Sample1
Sample2
```

*Execution:

For just 1 study use this snakemake command. Make sure your STAR_THREADS matches whats in -c. 
```
snakemake -j 25 -s STAR.py --cluster "sbatch -t 1:00:00 -c 30 -N 1"
```


For multiple studies, create 2 files:

* SRP: List of unique study accession IDs.
```
ERP126405
ERP127339
SRP293106
```
* list: 2 column file of study accession IDs and corresponding run accession IDs.
```
ERP124749       ERR4777044
ERP124749       ERR4777043
ERP126405       ERR5104751
ERP126405       ERR5104750
```
Then run STAR.py on all studies using this script. This will make it so each study gets its own combined SJ.out.tab file for the 2nd pass. 
```
cat SRP | while read i; do 
	cat list | grep "$i" | awk '{print $2}' > RAids.txt
	snakemake -j 300 -k -s STAR_SRA --cluster "sbatch -t 8:00:00 -c 30 -N 1 -p RM-shared"
	rm output/all.SJ.out.tab
done

```





## Infer Ancestry
Performs GATK best practices workflow for RNAseq short variant discovery. Intersects input varient data with varaint data from the 1000 Genomes Project to gather common ancestry informative loci. Performs PCA on variant data via PLINK and SVM model is implemented for ancestry inference. 

*In the config.yaml file, specify the number of PCs to be used for the machine learning step (Default 20) and chromosomes to be genotyped (Default Chr1-22).
*For SRA/fastq mode, simply continue to execution.


*For using ones own bam files:
	*Make sure to index bam with "samtools index -b 2ndPass.Aligned.sortedByCoord.out.bam"
	*ids.txt must contain directory names for individual sample bam(s). Below would be (Sample1,Sample2) 

OutputDir/Sample1/2ndPass.Aligned.sortedByCoord.out.bam
OutputDir/Sample2/2ndPass.Aligned.sortedByCoord.out.bam

```
ids.txt
```
Sample1
Sample2




Split RAids.txt so snakemake doesnt stall. 
```
split -l 100 RAids.txt

ls *xa* | cat > splits

cat splits | while read i; do
	cat $i > RAids.txt
	snakemake -j 300 -k -s InferAncestry.py --cluster "sbatch -t 02:00:00  -c 7 -p RM-shared"
done
```




