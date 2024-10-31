#!/bin/bash

#SBATCH --partition=short                      # Partition to run in
#SBATCH -c 1                                 # Requested cores
#SBATCH --time=0-01:00                    # Runtime in D-HH:MM format
#SBATCH --mem=20M                           # Requested Memory
#SBATCH -o %j.out                            # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e %j.err                            # File to which STDERR will be written, including job ID (%j)
#SBATCH --mail-type=ALL                      # ALL email notification type
#SBATCH --mail-user=<your_email_here>          # Email to which notifications will be sent

module load gcc/9.2.0 samtools/1.15.1


for name in bam/*_sorted.bam; do
	(basename ${name} _sorted.bam) >> logs/experimental_counts.log
	samtools view -c -F 388 ${name} >> logs/experimental_counts.log
done

for name in bam/spike_in/*_sorted.bam; do
	(basename ${name} _sorted.bam) >> logs/spikein_counts.log
	samtools view -c -F 388 ${name} >> logs/spikein_counts.log
done

