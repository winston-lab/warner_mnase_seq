#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  5 16:41:04 2024

@author: jlwarner
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

with open(r"logs/experimental_counts.log", 'r') as fp:
    lala = len(fp.readlines())
    print(lala)


# Read in .log files from samtools script output to generate lists of aligned read counts
libraries = []
scer_counts = []
spom_counts =[]
for i in range(0, lala, 2):
    file = open('logs/experimental_counts.log', 
                mode='r')
    content = file.readlines()
    libraries.append(content[i].strip('\n'))
    scer_counts.append(int(content[i+1].strip('\n')))
    file.close()
    file = open('logs/spikein_counts.log', 
                mode='r')
    content = file.readlines()
    spom_counts.append(int(content[i+1].strip('\n')))
    file.close()
    
libraries
scer_counts
spom_counts

# Read lists into data frame
aligned_reads = pd.DataFrame(
    {
     "library": libraries,
     "scer_counts": scer_counts,
     "spom_counts": spom_counts
     }
    )

aligned_reads
aligned_reads[['spom_counts', 'scer_counts']]

# Do some useful math and create new dataframe columns
aligned_reads["total_counts"] = aligned_reads["spom_counts"] + aligned_reads["scer_counts"]
aligned_reads["spom / scer"] = aligned_reads["spom_counts"] / aligned_reads["scer_counts"]
aligned_reads["proportion_spom"] = aligned_reads["spom_counts"] / aligned_reads["total_counts"]
aligned_reads["proportion_scer"] = aligned_reads["scer_counts"] / aligned_reads["total_counts"]
aligned_reads["1 / proportion_spom"] = 1 / aligned_reads["proportion_spom"]


# Check the math looks ok
# aligned_reads[["spom / scer", "proportion_spom", "1 / proportion_spom"]]


# Plot a stacked boxplot of the proportions of reads mapped to Spom & Scer genomes
aligned_reads[["library", "proportion_spom", "proportion_scer"]].plot(
    x = 'library',
    xlabel = 'proportion reads mapped',
    xlim = (0,1),
    kind = 'barh',
    stacked = True,
    figsize = (5,10),
    legend = False,
    title = 'Proportion of reads mapped to S. pombe \nor S. cerevisiae genomes',
    )
plt.savefig('logs/proportion_reads_mapped.png', dpi=300)

## Change x-axis limit to get a better look at Spom values
# aligned_reads[["library", "proportion_spom", "proportion_scer"]].plot(
#     x = 'library',
#     xlabel = 'proportion reads mapped',
#     xlim = (0,0.1),
#     kind = 'barh',
#     stacked = True,
#     figsize = (5,10),
#     legend = False,
#     title = 'Proportion of reads mapped to S. pombe \nor S. cerevisiae genomes, zoomed in',
#     )



### Testing some math
# aligned_reads["spom_counts"]*aligned_reads["1 / proportion_spom"]
# aligned_reads["scer_counts"]*aligned_reads["1 / proportion_spom"]
# aligned_reads["proportion_spom"]*aligned_reads["1 / proportion_spom"]
# (aligned_reads["spom_counts"]*aligned_reads["1 / proportion_spom"])/((aligned_reads["spom_counts"])+(aligned_reads["scer_counts"]))
# (aligned_reads["scer_counts"]*aligned_reads["1 / proportion_spom"])/((aligned_reads["spom_counts"])+(aligned_reads["scer_counts"]))
# (aligned_reads["scer_counts"]*aligned_reads["1 / proportion_spom"])/(aligned_reads["total_counts"])
# (aligned_reads["spom_counts"]*aligned_reads["1 / proportion_spom"])/((aligned_reads["spom_counts"]*aligned_reads["1 / proportion_spom"])+(aligned_reads["scer_counts"]*aligned_reads["1 / proportion_spom"]))



### Testing spike-in normalization
# aligned_reads["spom_counts"].max()
# aligned_reads["spom_counts"].max()/aligned_reads["spom_counts"]
aligned_reads["spike-in fold change"] = aligned_reads["spom_counts"].max()/aligned_reads["spom_counts"]
aligned_reads["spike-in fold change"]

for i in range(len(aligned_reads['spike-in fold change'])):
    aligned_reads.loc[i,'spike-in_norm_float'] = format(aligned_reads.loc[i, 'spike-in fold change'], '.16f')
# aligned_reads["spom_counts"]*aligned_reads["spike-in fold change"]
# aligned_reads["scer_counts"]*aligned_reads["spike-in fold change"]

aligned_reads[["library", "spike-in_norm_float"]].to_csv(path_or_buf='logs/normalization_table.csv', index=False, header=False)

