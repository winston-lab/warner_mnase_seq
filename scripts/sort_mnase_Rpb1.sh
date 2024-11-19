srun --pty -t 0-2:0:0 -p interactive /bin/bash

module load gcc/9.2.0 python/3.9.14 bedtools/2.30.0

bedtools intersect -a genome/annotations/Scer_transcripts_w_verifiedORFs-nonoverlapping_Rpb1sorted.bed -b genome/annotations/Scer_transcripts_w_verifiedORFs_nonoverlapping_WT-37C_plusonenuc.bed > genome/annotations/testA.bed

bedtools intersect -a genome/annotations/Scer_transcripts_w_verifiedORFs_nonoverlapping_WT-37C_plusonenuc.bed -b genome/annotations/Scer_transcripts_w_verifiedORFs-nonoverlapping_Rpb1sorted.bed > genome/annotations/testB.bed
