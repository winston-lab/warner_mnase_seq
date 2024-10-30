#!/bin/bash

#SBATCH --partition=short                      # Partition to run in
#SBATCH -c 1                                 # Requested cores
#SBATCH --time=0-00:05                    # Runtime in D-HH:MM format
#SBATCH --mem=500M                           # Requested Memory
#SBATCH -o %j.out                            # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e %j.err                            # File to which STDERR will be written, including job ID (%j)
#SBATCH --mail-type=ALL                      # ALL email notification type
#SBATCH --mail-user=james_warner@hms.harvard.edu          # Email to which notifications will be sent

#Use older versions to avoid incompatibility with python 'dpp_inches'
module load gcc/6.2.0 python/2.7.12 deeptools/3.0.2	
	

plotCorrelation -in correlation/scores_per_bin.gz \
	-c spearman \
	-p heatmap \
	-o correlation/spearman.png \
	--skipZeros \
	--plotNumbers \
	--outFileCorMatrix correlation/spearman.tab

plotCorrelation -in correlation/scores_per_bin.gz \
	-c pearson \
	-p heatmap \
	-o correlation/pearson.png \
	--skipZeros \
	--plotNumbers \
	--outFileCorMatrix correlation/pearson.tab

