#!/bin/zsh

conda activate bioinf

###############
## PART 2
###############

# create a project directory on your computer in your GitHub directory. Name this directory 'seqaln' for sequence alignment.
# and set the 'PROJ_DIR' variable to the path to our working directory

mkdir seqaln
export PROJ_DIR=/Users/ojohnson/Documents/GitHub/Bioinformatics_Spring2024/seqaln

###############
## Genome setup
###############

#Get genome files
# first, create the reference genome directory and download the genome files into it

mkdir -p $PROJ_DIR/genome && cd $PROJ_DIR/genome

# next, use the wget command to 'fetch' the genome files from the ncbi website. These will be saved to your current working directory (in this case the 'genome' folder). We'll use this one: 
# https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_900604845.1/

# download it using the wget command, which dowloads stuff off the internet
wget https://www.ncbi.nlm.nih.gov/assembly/GCA_900604845.1

#unzip the DNA sequence & index
unzip ncbi_dataset.zip 

# create a symbolic link with a simpler name for the genome
ln -s ncbi_dataset/data/GCA_900604845.1/GCA_900604845.1_TTHNAR1_genomic.fna Thermus_thermophilus_TTHNAR1.fa

# How many sequences are in the genome assembly for this bacterium?
grep ">" file.fasta | wc -l


#samtools index genome
# https://www.biostars.org/p/212594/
# Indexing a genome can be explained similar to indexing a book. 
# If you want to know on which page a certain word appears or a chapter begins, 
# it is much more efficient/faster to look it up in a pre-built index than going 
# through every page of the book until you found it. Same goes for alignments. 
# Indices allow the aligner to narrow down the potential origin of a query sequence 
# within the genome, saving both time and memory.

samtools faidx Thermus_thermophilus_TTHNAR1.fa

# picard dictionary for genome

picard CreateSequenceDictionary \
    REFERENCE=Thermus_thermophilus_TTHNAR1.fa \
    OUTPUT=Thermus_thermophilus_TTHNAR1.dict

# build bowtie2 genome index

bowtie2-build Thermus_thermophilus_TTHNAR1.fa Thermus_thermophilus_TTHNAR1

#####################
## Get sequence reads
#####################

# https://trace.ncbi.nlm.nih.gov/Traces/sra/?run=SRR5324768
# from "Recovery of nearly 8,000 metagenome-assembled genomes substantially expands the tree of life."
# Nat Microbiol. 2017 Nov;2(11):1533-1542. doi: 10.1038/s41564-017-0012-7. Epub 2017 Sep 11.

# get reads from SRA, put into project directory

cd $PROJ_DIR

mkdir fastq && cd fastq

# assign sample number to a variable called SRR
export SRR=SRR5324768

# download forward and reverse reads
wget https://github.com/marctollis/INF515-Comparative-Genomics_fall22/raw/main/labs/lab_1/${SRR}_pass_1.fastq.gz
wget https://github.com/marctollis/INF515-Comparative-Genomics_fall22/raw/main/labs/lab_1/${SRR}_pass_2.fastq.gz

# if you want to rename them to something shorter:
# mv ${SRR}_pass_1.fastq.gz ${SRR}_1.fastq.gz
# mv ${SRR}_pass_2.fastq.gz ${SRR}_2.fastq.gz

# look at the first few lines of each file (a few options)
gzcat ${SRR}_pass_1.fastq.gz | head
gzcat ${SRR}_pass_2.fastq.gz | head 

# gzip -cd ${SRR}_pass_1.fastq.gz | head
# gzip -cd ${SRR}_pass_2.fastq.gz | head

# less -N ${SRR}_pass_1.fastq.gz
# less -N ${SRR}_pass_2.fastq.gz


#####################
## Alignment Time!!!!
#####################

cd $PROJ_DIR

mkdir -p alignment

# http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml
# http://www.htslib.org/doc/samtools.html

# this step takes about 5 minutes 

bowtie2 -x genome/Thermus_thermophilus_TTHNAR1 \
        -1 fastq/${SRR}_pass_1.fastq.gz \
        -2 fastq/${SRR}_pass_2.fastq.gz --sensitive-local \
        --rg-id ${SRR} --rg SM:${SRR} --rg PL:ILLUMINA \
    | samtools view -hb - | samtools sort -l 5 -o alignment/${SRR}.bam

samtools index alignment/${SRR}.bam

# text view alignment with 
samtools tview alignment/SRR5324768.bam genome/Thermus_thermophilus_TTHNAR1.fa

# type ? for help, q to quit
# pileup format - is a text-based format for summarizing the base calls of aligned reads to a reference sequence.

##########################
## Variant Calls with GATK
##########################
# https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php
#This step may take 5-20 minutes depending on your hardware

cd $PROJ_DIR

mkdir -p variants

# we will give this job a little more RAM than the default since gatk gets hungry. 
# this part will typically take ~10-15 minutes.

gatk --java-options "-Xmx8g" HaplotypeCaller  \
   --reference genome/Thermus_thermophilus_TTHNAR1.fa \
   --sample-ploidy 1 \
   --input alignment/${SRR}.bam \
   --output variants/${SRR}.vcf

# the output is in the variant call format (vcf)
# https://samtools.github.io/hts-specs/VCFv4.2.pdf
# VCF is a text file format (most likely stored in a compressed manner). 
# It contains meta-information lines, a header line, and then data lines 
# each containing information about a position in the genome. 
# The format also has the ability
# to contain genotype information on samples for each position.

# https://gatk.broadinstitute.org/hc/en-us/articles/360035531692-VCF-Variant-Call-Format

less variants/${SRR}.vcf


conda deactivate

