
# RNA-Seq Inferred Ancestry (RIA)
RIA is a method for infering the super-population(Africa, Europe, South Asia, East Asia, and America) from Human RNA-seq data.



![alt text](https://github.com/jahaltom/RNA-Seq-Ancestry-Inference/blob/main/FlowChart.png?raw=true)

## Prerequisites

Using Conda 4.10.3, create the conda enviroment and activate:
```
conda env create -f environment.yml
conda activate Ancestry
```


## Data Preparation

**1000 Genomes Project:**
The snakemake script "Prepare_1KGP" downloads chr(1-22) level VCF files from 1000 Genomes Project phase 3 on GRCh38 (https://www.internationalgenome.org/data-portal/data-collection/grch38, https://doi.org/10.12688/wellcomeopenres.15126.2) and filters out INDELs along with indexing the resulting VCF. It also creates the interval lists needed fot the analysis. 
```
snakemake -j 22 -s Prepare_1KGP --cluster "sbatch -t 00:30:00 -c 4 -N 1"
```

**GRCh38 Reference Genome**
The bash script "Prepare_Reference_Genome" will download the Human genome GRCh38 fasta(GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz) and the corresponding gtf, and will create a seqence dictionary and index file for the fasta. It also creates a STAR index.
```
sbatch Prepare_Reference_Genome
```

## Raw data retrieval from SRA, QC, and STAR 2-Pass

The snakemake script "STAR_SRA" takes in a list of run accession IDs "RAids.txt" and fetches the raw fastq files from SRA and then uses Trimgalore for QC. The reads are then ran through STAR 2-Pass mode for enhanced novel SJ detection. The SJ.out.tab file for the 2nd pass is made by combining all SJ.out.tab files from the first pass and removing SJ's that are supported by 2 or less unique mappers. 

For just 1 study, create a list of the corresponding run accession IDs "RAids.txt" and run
```
snakemake -j 50 -k -s STAR_SRA --cluster "sbatch -t 8:00:00 -c 30 -N 1"
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
Then run STAR_SRA on all studies using this script. This will make it so each study gets its own combined SJ.out.tab file for the 2nd pass. 
```
cat SRP | while read i; do 
	cat list | grep "$i" | awk '{print $2}' > RAids.txt
	snakemake -j 50 -k -s STAR_SRA --cluster "sbatch -t 8:00:00 -c 30 -N 1 -p biocrunch"
	rm output/all.SJ.out.tab
done

```








```
export OMP_NUM_THREADS=3
export GIT_PYTHON_REFRESH=quiet 
```


