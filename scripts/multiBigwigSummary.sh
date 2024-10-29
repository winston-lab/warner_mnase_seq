#!/bin/bash

#SBATCH --partition=short                      # Partition to run in
#SBATCH -c 4                                 # Requested cores
#SBATCH --time=0-00:20                    # Runtime in D-HH:MM format
#SBATCH --mem=500M                           # Requested Memory
#SBATCH -o %j.out                            # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e %j.err                            # File to which STDERR will be written, including job ID (%j)
#SBATCH --mail-type=ALL                      # ALL email notification type
#SBATCH --mail-user=<YOUR_EMAIL_HERE>          # Email to which notifications will be sent


module load gcc/9.2.0 python/3.9.14 deeptools/3.5.0

multiBigwigSummary bins -b deeptools/*.bw \
	-o correlation/scores_per_bin.gz \
	--smartLabels \
	-bs 75 \
	-p max \
	--chromosomesToSkip chrM \


# Add following argument if you want a table of raw scores:
# 	--outRawCounts correlation/scores_per_bin.tab