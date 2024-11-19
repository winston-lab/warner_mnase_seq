#!/bin/bash

#SBATCH --partition=short                      # Partition to run in
#SBATCH -c 4                                 # Requested cores
#SBATCH --time=0-00:10                   # Runtime in D-HH:MM format
#SBATCH --mem=300M                           # Requested Memory
#SBATCH -o %j.out                            # File to which STDOUT will be written, including job ID (%j)
#SBATCH -e %j.err                            # File to which STDERR will be written, including job ID (%j)
#SBATCH --mail-type=ALL                      # ALL email notification type
#SBATCH --mail-user=<your_email_here>          # Email to which notifications will be sent


module load gcc/9.2.0 python/3.9.14 deeptools/3.5.0


computeMatrix reference-point -S deeptools/*${1%}*.bw -R genome/annotations/Scer_transcripts_w_verifiedORFs_nonoverlapping_WT-37C_plusonenuc_Rpb1sorted.bed -o deeptools/${1%}_plusone_sorted.gz \
	-b 500 \
	-a 1500 \
	-bs 10 \
	-p max \
	--averageTypeBins mean \
	--outFileNameMatrix deeptools/tab/${1%}_plusone_sorted.tab
	