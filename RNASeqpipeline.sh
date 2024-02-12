#!/bin/bash

SECONDS=0

# change the working directory
cd transcriptome

# STEP 1: Run fastqc
fastqc -t 4 raw_data/*.fastq -o raw_data/

# run trimmomatic to trim reads with poor quality and also --fastqc option will run fastqc automatically
trim_galore --fastqc --core 2 raw_data/SRR453566.fastq --phred33 -o filtered_data
echo "Trim_galore finished running!"

##running and showing demo with full file will take time that's why making subset of filtered fastq file
##now make a subset of trim fastq file (fetch 0.5 million reads from filtered fastq file) so that you can show the mapping.
head -500000 filtered_data/SRR453566_trimmed.fq > filtered_data/SRR453566_trimmed_0.5M.fq

# mkdir HISAT2
# get the genome
#wget -c http://ftp.ensembl.org/pub/release-107/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz

## after downloading above fa.gz
# gunzip Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz
# mv Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa genome.fa

## get the gtf annotation file
# wget –c http://ftp.ensembl.org/pub/release-107/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.107.gtf.gz
# gunzip Saccharomyces_cerevisiae.R64-1-1.107.gtf.gz
# mv Saccharomyces_cerevisiae.R64-1-1.107.gtf genome.gtf

##make a file which contain splice site

hisat2_extract_splice_sites.py genome.gtf > splice_sites.txt

##create a genome index file
hisat2-build genome.fa genome_index

# STEP 2: Run HISAT2

# create a folder to store alignment
mkdir alignment

#Run HISAT2
hisat2 -x genome_index --known-splicesite-infile splice_sites.txt -p 2 -U filtered_data/SRR453566_trimmed_0.5M.fq -S alignment/SRR453566_trimmed_0.5M.out.sam

##convert sam alignment file into bam file
samtools view -bS alignment/SRR453566_trimmed_0.5M.out.sam | samtools sort > alignment/SRR453566_trimmed_0.5M.out_sorted.bam

##hisat2 -q --rna-strandness R -x genome -U filtered_data/demo_trimmed.fastq | samtools sort -o /demo_trimmed.bam
echo "HISAT2 finished running!"

# STEP 3: Run featureCounts - Quantification: this will be used in differential gene expression analysis
mkdir Quantification
subread-2.0.6-Windows-x86_64/bin/featureCounts.exe -T 2 -a genome.gtf -o Quantification/SRR453566.featureCounts.txt \
alignment/SRR453566_trimmed_0.5M.out_sorted.bam

echo "featureCounts finished running!"

#STEP 4: Estimate transcript/gene expression using stringtie this could be used in heatmap
mkdir expression
stringtie/stringtie alignment/SRR453566_trimmed_0.5M.out_sorted.bam -o expression/transcripts.gtf -e -G genome.gtf -A expression/gene_abund.tab
duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."

