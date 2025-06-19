#!/bin/bash

#SBATCH --partition=short                      # Partition to run in
#SBATCH -c 4                                 # Requested cores
#SBATCH --time=0-03:00                    # Runtime in D-HH:MM format
#SBATCH --mem=2GB                           # Requested Memory
#SBATCH -o %j.out                            # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e %j.err                            # File to which STDERR will be written, including job ID (%j)
#SBATCH --mail-type=ALL                      # ALL email notification type
#SBATCH --mail-user=<YOUR_EMAIL_HERE>          # Email to which notifications will be sent

module load gcc/6.2.0 bowtie2/2.2.9 samtools/1.3.1

dir=/n/groups/winston/jw362/mnase_seq_new/fastq
base=$(basename ${1%} _R1_001.fastq.gz)

bowtie2 --no-unal \
	--no-mixed \
	--no-discordant \
	-q \
	--sensitive \
	-p 4 \
	-x genome/bowtie2_index/S_cerevisiae.R64-2-1 \
	-1 $dir/${base}_R1_001.fastq.gz -2 $dir/${base}_R2_001.fastq.gz \
	2> logs/${base}_bowtie2.txt \
	| samtools view - -bS -q 30 > bam/${base}_unsorted.bam

samtools sort bam/${base}_unsorted.bam -o bam/${base}_sorted.bam

samtools index bam/${base}_sorted.bam

