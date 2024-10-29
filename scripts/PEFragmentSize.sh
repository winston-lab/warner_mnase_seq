#!/bin/bash

#SBATCH --partition=short                      # Partition to run in
#SBATCH -c 1                                 # Requested cores
#SBATCH --time=0-00:05                    # Runtime in D-HH:MM format
#SBATCH --mem=150M                           # Requested Memory
#SBATCH -o %j.out                            # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e %j.err                            # File to which STDERR will be written, including job ID (%j)
#SBATCH --mail-type=ALL                      # ALL email notification type
#SBATCH --mail-user=james_warner@hms.harvard.edu          # Email to which notifications will be sent

module load gcc/9.2.0 python/3.9.14 deeptools/3.5.0

bamPEFragmentSize -b bam/*_sorted.bam -o fragment_sizes/PEF_histogram.png --table fragment_sizes/PEF_table.csv

bamPEFragmentSize -b bam/spike_in/*_sorted.bam -o fragment_sizes/spikein_PEF_histogram.png --table fragment_sizes/spikein_PEF_table.csv
