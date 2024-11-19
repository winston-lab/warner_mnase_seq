#!/bin/bash

#SBATCH --partition=short                      # Partition to run in
#SBATCH -c 1                                 # Requested cores
#SBATCH --time=0-00:05                    # Runtime in D-HH:MM format
#SBATCH --mem=750M                           # Requested Memory
#SBATCH -o %j.out                            # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e %j.err                            # File to which STDERR will be written, including job ID (%j)
#SBATCH --mail-type=ALL                      # ALL email notification type
#SBATCH --mail-user=<your_email_here>          # Email to which notifications will be sent

#Use older versions to avoid incompatibility with python 'dpp_inches'
module load gcc/6.2.0 python/2.7.12 deeptools/3.0.2


plotHeatmap -m deeptools/*${1%}*_plusone_sorted.gz -o plots/${1%}_plusone_sorted_heatmap.png \
	--dpi 300 \
	--refPointLabel "+1 dyad" \
        --heatmapWidth 4 \
	--heatmapHeight 5 \
	--whatToShow "heatmap and colorbar" \
	--sortRegions keep \
	--xAxisLabel "position" \
	--outFileNameMatrix plots/${1%}_plusone_sorted_heatmap.tab \

plotHeatmap -m deeptools/*${1%}*_plusone_sorted.gz -o plots/${1%}_plusone_sorted_heatmap_limit.png \
	--dpi 300 \
	--refPointLabel "+1 dyad" \
        --heatmapWidth 4 \
	--heatmapHeight 5 \
	--whatToShow "heatmap and colorbar" \
	--sortRegions keep \
	--xAxisLabel "position" \
	-max 30 \

# If you would like to sort using a particular sample, add the following line of code:
#	--sortUsingSamples <SAMPLE NUMBER> \
