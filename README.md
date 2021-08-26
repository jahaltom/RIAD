**Prerequisites**
Create the conda enviroment and activate:
```
conda env create -f environment.yml
conda activate Ancestry
```
Activate GATK4"
```
export PATH="/work/LAS/xgu-lab/tools/gatk-4.2.0.0/:$PATH"
```

**Data Preparation**
This snakemake script creates a directory called "data" and then takes the 1000 Genomes Project Phase3 on GRCh38 VCF files(Individual Chr) from the web without needed to download and filters out INDELs creating a filtered VCF along with indexing the resulting VCF. It also creates an interval list for each Chr to be used with GATK4. 
* Prepare 1000G data
