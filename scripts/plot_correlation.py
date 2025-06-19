#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 30 16:23:24 2025

@author: jlwarner
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy

renamer = {
    "'425_414_I_rep3_S26_si'": 'EV_POB3_IAA_rep3',
    "'425_95_I_rep1_S6_si'": 'Δ2-238_POB3_IAA_rep1',
    "'425_95_I_rep2_S18_si'": 'Δ2-238_POB3_IAA_rep2',
    "'425_95_I_rep3_S30_si'": 'Δ2-238_POB3_IAA_rep3',
    "'425_414_I_rep1_S2_si'": 'EV_POB3_IAA_rep1',
    "'425_414_I_rep2_S14_si'": 'EV_POB3_IAA_rep2',
    "'553_414_I_rep1_S8_si'": 'EV_E154K_IAA_rep1',
    "'553_414_I_rep2_S20_si'": 'EV_E154K_IAA_rep2',
    "'553_414_I_rep3_S32_si'": 'EV_E154K_IAA_rep3',
    "'553_95_I_rep3_S36_si'": 'Δ2-238_E154K_IAA_rep3',
    "'553_95_I_rep1_S12_si'": 'Δ2-238_E154K_IAA_rep1',
    "'553_95_I_rep2_S24_si'": 'Δ2-238_E154K_IAA_rep2',
    "'425_93_I_rep1_S4_si'": 'SPT6_POB3_IAA_rep1',
    "'425_414_D_rep3_S25_si'": 'EV_POB3_DMSO_rep3',
    "'425_93_I_rep2_S16_si'": 'SPT6_POB3_IAA_rep2',
    "'553_93_I_rep2_S22_si'": 'SPT6_E154K_IAA_rep2',
    "'553_93_I_rep1_S10_si'": 'SPT6_E154K_IAA_rep1',
    "'553_93_I_rep3_S34_si'": 'SPT6_E154K_IAA_rep3',
    "'425_95_D_rep1_S5_si'": 'Δ2-238_POB3_DMSO_rep1',
    "'425_93_I_rep3_S28_si'": 'SPT6_POB3_IAA_rep3',
    "'425_95_D_rep2_S17_si'": 'Δ2-238_POB3_DMSO_rep2',
    "'425_93_D_rep1_S3_si'": 'SPT6_POB3_DMSO_rep1',
    "'425_93_D_rep2_S15_si'": 'SPT6_POB3_DMSO_rep2',
    "'553_95_D_rep3_S35_si'": 'Δ2-238_E154K_DMSO_rep3',
    "'553_414_D_rep2_S19_si'": 'EV_E154K_DMSO_rep2',
    "'553_95_D_rep2_S23_si'": 'Δ2-238_E154K_DMSO_rep2',
    "'553_93_D_rep2_S21_si'": 'SPT6_E154K_DMSO_rep2',
    "'553_93_D_rep3_S33_si'": 'SPT6_E154K_DMSO_rep3',
    "'553_93_D_rep1_S9_si'": 'SPT6_E154K_DMSO_rep1',
    "'553_95_D_rep1_S11_si'": 'Δ2-238_E154K_DMSO_rep1',
    "'553_414_D_rep1_S7_si'": 'EV_E154K_DMSO_rep1',
    "'553_414_D_rep3_S31_si'": 'EV_E154K_DMSO_rep3',
    "'425_414_D_rep1_S1_si'": 'EV_POB3_DMSO_rep1',
    "'425_414_D_rep2_S13_si'": 'EV_POB3_DMSO_rep2',
    "'425_93_D_rep3_S27_si'": 'SPT6_POB3_DMSO_rep3',
    "'425_95_D_rep3_S29_si'": 'Δ2-238_POB3_DMSO_rep3',
    }

# plotting correlation
pearson_corr = pd.read_csv('/Users/jlwarner/Desktop/mnase_seq_new/correlation/tab/mnase_pearson.tab', 
                 sep='\t', header=1)
pearson_corr.rename(columns=renamer,inplace=True)
pearson_corr.rename(index=renamer,inplace=True)
sns.set_style('ticks')
plt.figure(dpi=300)
g = sns.clustermap(pearson_corr, 
               xticklabels=True, yticklabels=True, figsize=(10,10),
               annot=False, 
#               annot_kws={'fontweight':'bold'},
               )
g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(), fontweight='bold', fontsize=10)
g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontweight='bold', fontsize=10)
plt.savefig('plots/correlation/pearson_corr.png', dpi=300)
plt.show();


