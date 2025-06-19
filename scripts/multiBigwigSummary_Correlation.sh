#!/bin/bash

#SBATCH --partition=short                      # Partition to run in
#SBATCH -c 4                                 # Requested cores
#SBATCH --time=0-00:10                    # Runtime in D-HH:MM format
#SBATCH --mem=400M                           # Requested Memory
#SBATCH -o %j.out                            # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e %j.err                            # File to which STDERR will be written, including job ID (%j)
#SBATCH --mail-type=ALL                      # ALL email notification type
#SBATCH --mail-user=<YOUR_EMAIL_HERE>          # Email to which notifications will be sent


module load gcc/9.2.0 python/3.10.11

source env/deeptools/bin/activate

multiBigwigSummary bins -b deeptools/si/*_si.bw \
	-o correlation/gz/mnase_scores_per_bin.gz \
	--smartLabels \
	-bs 75 \
	-p max \
	--chromosomesToSkip chrM \

plotCorrelation -in correlation/gz/mnase_scores_per_bin.gz \
	-c spearman \
	-p heatmap \
	-o correlation/plots/mnase_spearman.png \
	--skipZeros \
	--plotNumbers \
	--outFileCorMatrix correlation/tab/mnase_spearman.tab

plotCorrelation -in correlation/gz/mnase_scores_per_bin.gz \
	-c pearson \
	-p heatmap \
	-o correlation/plots/mnase_pearson.png \
	--skipZeros \
	--plotNumbers \
	--outFileCorMatrix correlation/tab/mnase_pearson.tab

deactivate


# Add following argument if you want a table of raw scores:
# 	--outRawCounts correlation/scores_per_bin.tab