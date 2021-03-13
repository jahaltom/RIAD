#download human transcriptome
wget -q ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/gencode.v36.transcripts.fa.gz
#human genome
wget -q ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
#download covid seq
wget -q ftp://ftp.ensemblgenomes.org/pub/viruses/fasta/sars_cov_2/cdna/Sars_cov_2.ASM985889v3.cdna.all.fa.gz
#Download GTF
wget http://ftp.ensembl.org/pub/release-103/gtf/homo_sapiens/Homo_sapiens.GRCh38.103.chr_patch_hapl_scaff.gtf.gz
gunzip -f *.gz


##Extraxt scaffolds
grep ">" GCA_000001405.15_GRCh38_no_alt_analysis_set.fna | awk '{print $1}' | sed 's/>//g' > HumanScaffolds.txt
#combine transcriptome and decoy fasta files
cat gencode.v36.transcripts.fa Sars_cov_2.ASM985889v3.cdna.all.fa GCA_000001405.15_GRCh38_no_alt_analysis_set.fna > human_tr_gen_decoy.fasta

#create salmon index
time salmon index -t human_tr_gen_decoy.fasta -d HumanScaffolds.txt -p 15 -i salmon_index






















####Installation

#Load conda
module load  miniconda3/4.3.30-qdauveb

conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge


#Create a new conda enviroment and activate
conda create -n pyrpipe python=3.8
source activate pyrpipe
#Install tools 
conda install -c bioconda pyrpipe star=2.7.7a sra-tools=2.10.9 stringtie=2.1.4 trim-galore=0.6.6 orfipy=0.0.3 salmon=1.4.0

#To create a yaml file containing information about the conda environment, run the following command
conda env export | grep -v "^prefix: " > environment.yml

#To recreate the conda environment in the environment.yml, use

conda env create -f environment.yml

module load snakemake/3.11.2-py3-wait26e



from pyrpipe import sra,qc,mapping,assembly
#define some vaiables
run_id='SRR976159'
working_dir='example_output'
gen='Arabidopsis_thaliana.TAIR10.dna.toplevel.fa'
ann='Arabidopsis_thaliana.TAIR10.46.gtf'
star_index='star_index/athaliana'
#initialize objects
#creates a star object to use with threads
star=mapping.Star(index=star_index,genome=gen,threads=4)
#use trim_galore for trimming
trim_galore=qc.Trimgalore()
#Stringtie for assembly
stringtie=assembly.Stringtie(guide=ann)
#create SRA object; this will download fastq if doesnt exist
srr_object=sra.SRA(run_id,directory=working_dir)
#create a pipeline using the objects
srr_object.trim(trim_galore).align(star).assemble(stringtie)

#The assembled transcripts are in srr_object.gtf
print('Final result',srr_object.gtf)


