#!/bin/bash

#SBATCH --partition=short                      # Partition to run in
#SBATCH -c 4                                 # Requested cores
#SBATCH --time=0-02:00                    # Runtime in D-HH:MM format
#SBATCH --mem=2GB                           # Requested Memory
#SBATCH -o %j.out                            # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e %j.err                            # File to which STDERR will be written, including job ID (%j)
#SBATCH --mail-type=ALL                      # ALL email notification type
#SBATCH --mail-user=james_warner@hms.harvard.edu          # Email to which notifications will be sent

module load gcc/6.2.0 bowtie2/2.2.9 samtools/1.3.1 java/jdk-1.8u112 qualimap/2.2.1 bedtools/2.27.1

base=$(basename ${1%} _R1_001.fastq.gz)

bowtie2 --no-unal \
	-q \
	-k 2 \
	-p 4 \
	-x genome/bowtie2_index/S_cerevisiae.R64-2-1 \
	-1 fastq/${base}_R1_001.fastq.gz -2 fastq/${base}_R2_001.fastq.gz \
	2> logs/${base}_bowtie2.txt \
	| samtools view - -bS -q 30 > bam/${base}_unsorted.bam

samtools sort bam/${base}_unsorted.bam -o bam/${base}_sorted.bam

samtools index bam/${base}_sorted.bam
