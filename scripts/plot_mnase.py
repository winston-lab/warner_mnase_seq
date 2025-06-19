#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  4 21:21:22 2025

@author: jlwarner
"""


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scipy




#Sorted by Rpb1 enrichment
positions = np.arange(-500,1500,10)

spt6_pob3_iaa_averaged = pd.read_csv('/Users/jlwarner/Desktop/mnase_seq_new/transfer/tab/425_93_I_averaged_plusone_sorted.tab', 
            sep='\t', header=3, names=positions)
ev_pob3_iaa_averaged = pd.read_csv('/Users/jlwarner/Desktop/mnase_seq_new/transfer/tab/425_414_I_averaged_plusone_sorted.tab', 
            sep='\t', header=3, names=positions)
dntd_pob3_iaa_averaged = pd.read_csv('/Users/jlwarner/Desktop/mnase_seq_new/transfer/tab/425_95_I_averaged_plusone_sorted.tab', 
            sep='\t', header=3, names=positions)
spt6_e154k_iaa_averaged = pd.read_csv('/Users/jlwarner/Desktop/mnase_seq_new/transfer/tab/553_93_I_averaged_plusone_sorted.tab', 
            sep='\t', header=3, names=positions)
ev_e154k_iaa_averaged = pd.read_csv('/Users/jlwarner/Desktop/mnase_seq_new/transfer/tab/553_414_I_averaged_plusone_sorted.tab', 
            sep='\t', header=3, names=positions)
dntd_e154k_iaa_averaged = pd.read_csv('/Users/jlwarner/Desktop/mnase_seq_new/transfer/tab/553_95_I_averaged_plusone_sorted.tab', 
            sep='\t', header=3, names=positions)

spt6_pob3_dmso_averaged = pd.read_csv('/Users/jlwarner/Desktop/mnase_seq_new/transfer/tab/425_93_D_averaged_plusone_sorted.tab', 
            sep='\t', header=3, names=positions)
ev_pob3_dmso_averaged = pd.read_csv('/Users/jlwarner/Desktop/mnase_seq_new/transfer/tab/425_414_D_averaged_plusone_sorted.tab', 
            sep='\t', header=3, names=positions)
dntd_pob3_dmso_averaged = pd.read_csv('/Users/jlwarner/Desktop/mnase_seq_new/transfer/tab/425_95_D_averaged_plusone_sorted.tab', 
            sep='\t', header=3, names=positions)
spt6_e154k_dmso_averaged = pd.read_csv('/Users/jlwarner/Desktop/mnase_seq_new/transfer/tab/553_93_D_averaged_plusone_sorted.tab', 
            sep='\t', header=3, names=positions)
ev_e154k_dmso_averaged = pd.read_csv('/Users/jlwarner/Desktop/mnase_seq_new/transfer/tab/553_414_D_averaged_plusone_sorted.tab', 
            sep='\t', header=3, names=positions)
dntd_e154k_dmso_averaged = pd.read_csv('/Users/jlwarner/Desktop/mnase_seq_new/transfer/tab/553_95_D_averaged_plusone_sorted.tab', 
            sep='\t', header=3, names=positions)




ticks = [-500,0,500,1000,1499]
labels = ['-0.5 kb', '+1 dyad', '+0.5 kb', '+1 kb','+1.5 kb']

#metagenes
sns.set_style('ticks', {'axes.facecolor':'white'})

## Figure 5 and S5
plt.figure(dpi=300)
plt.plot(positions, spt6_pob3_iaa_averaged.mean(), label = 'Spt6', color = 'black', lw=1.5)
plt.plot(positions, ev_pob3_iaa_averaged.mean(), label = 'empty vector', color = 'mediumblue', lw=1.5)
plt.plot(positions, dntd_pob3_iaa_averaged.mean(), label = 'Spt6Δ2-238', color = 'red', lw=1.5)
plt.xlabel('position', fontsize=14)
plt.xlim(left=-250, right=1250)
plt.xticks(ticks=ticks, labels=labels, fontsize=12)
plt.ylabel('dyad occupancy', fontsize=14)
plt.ylim(0, 30)
plt.yticks(ticks=np.array((0, 10,20,30)), fontsize=12)
plt.legend(loc='upper right', fontsize=12)
#plt.title('POB3 IAA')
plt.tight_layout()
plt.savefig('plots/metagenes/fig5A.svg')
plt.show();

plt.figure(dpi=300)
plt.plot(positions, spt6_e154k_iaa_averaged.mean(), label = 'Spt6\nPob3-E154K', color = 'black', lw=1.5)
plt.plot(positions, dntd_pob3_iaa_averaged.mean(), label = 'Spt6Δ2-238', color = 'red', lw=1.5)
plt.plot(positions, dntd_e154k_iaa_averaged.mean(), label = 'Spt6Δ2-238\nPob3-E154K', color = 'mediumorchid', lw=1.5)
plt.xlabel('position', fontsize=14)
plt.xlim(left=-250, right=1250)
plt.xticks(ticks=ticks, labels=labels, fontsize=12)
plt.ylabel('dyad occupancy', fontsize=14)
plt.ylim(0, 40)
plt.yticks(ticks=np.array((0, 10,20,30,40)), fontsize=12)
plt.legend(loc='upper right', fontsize=12)
plt.tight_layout()
plt.savefig('plots/metagenes/fig5D.svg')
plt.show();

plt.figure(dpi=300)
plt.plot(positions, spt6_e154k_iaa_averaged.mean(), label = 'Spt6\nPob3-E154K', color = 'black', lw=1.5)
plt.plot(positions, ev_pob3_iaa_averaged.mean(), label = 'empty vector', color = 'mediumblue', lw=1.5)
plt.plot(positions, ev_e154k_iaa_averaged.mean(), label = 'empty vector\nPob3-E514K', color = 'mediumaquamarine', lw=1.5)
plt.xlabel('position', fontsize=14)
plt.xlim(left=-250, right=1250)
plt.xticks(ticks=ticks, labels=labels, fontsize=12)
plt.ylabel('dyad occupancy', fontsize=14)
plt.ylim(0, 40)
plt.yticks(ticks=np.array((0, 10,20,30,40)), fontsize=12)
plt.legend(loc='upper right', fontsize=12)
#plt.title('POB3 IAA')
plt.tight_layout()
plt.savefig('plots/metagenes/fig5E.svg')
plt.show();


plt.figure(dpi=300)
plt.plot(positions, spt6_pob3_dmso_averaged.mean(), label = 'Spt6 DMSO', color = 'black', lw=1.5)
plt.plot(positions, spt6_pob3_iaa_averaged.mean(), label = 'Spt6 IAA', color = 'crimson', lw=1.5)
plt.xlabel('position', fontsize=14)
plt.xlim(left=-250, right=1250)
plt.xticks(ticks=ticks, labels=labels, fontsize=12)
plt.ylabel('dyad occupancy', fontsize=14)
plt.ylim(0, 30)
plt.yticks(ticks=np.array((0, 10,20,30)), fontsize=12)
plt.legend(loc='upper right', fontsize=12)
#plt.title('POB3 IAA')
plt.tight_layout()
plt.savefig('plots/metagenes/figS5C.svg')
plt.show();

plt.figure(dpi=300)
plt.plot(positions, spt6_pob3_iaa_averaged.mean(), label = 'Spt6', color = 'black', lw=1.5)
plt.plot(positions, spt6_e154k_iaa_averaged.mean(), label = 'Spt6\nPob3-E154K', color = 'orangered', lw=1.5)
plt.xlabel('position', fontsize=14)
plt.xlim(left=-250, right=1250)
plt.xticks(ticks=ticks, labels=labels, fontsize=12)
plt.ylabel('dyad occupancy', fontsize=14)
plt.ylim(0, 40)
plt.yticks(ticks=np.array((0, 10,20,30,40)), fontsize=12)
plt.legend(loc='upper right', fontsize=12)
#plt.title('POB3 IAA')
plt.tight_layout()
plt.savefig('plots/metagenes/figS5D.svg')
plt.show();


###Other metagene plots
plt.figure(dpi=300)
plt.plot(positions, spt6_e154k_iaa_averaged.mean(), label = 'SPT6 IAA', color = 'black', lw=2)
plt.plot(positions, dntd_e154k_iaa_averaged.mean(), label = 'ΔNTD IAA', color = 'red', lw=2)
plt.plot(positions, ev_e154k_iaa_averaged.mean(), label = 'EV IAA', color = 'mediumblue', lw=2)
plt.xlabel('position')
plt.xlim(left=-250, right=1250)
plt.xticks(ticks=ticks, labels=labels)
plt.ylabel('spike-in normalized dyad occupancy')
plt.ylim(0, 40)
plt.yticks(ticks=np.array((0, 10,20,30,40)))
plt.legend(loc='upper right')
plt.title('E154K')
plt.tight_layout()
plt.savefig('plots/metagenes/E154K_IAA.png', dpi=300)
plt.show();


plt.figure(dpi=300)
plt.plot(positions, spt6_pob3_dmso_averaged.mean(), label = 'SPT6 DMSO', color = 'black', lw=2)
plt.plot(positions, dntd_pob3_dmso_averaged.mean(), label = 'ΔNTD DMSO', color = 'red', lw=2)
plt.plot(positions, ev_pob3_dmso_averaged.mean(), label = 'EV DMSO', color = 'darkorange', lw=2)
plt.xlabel('position')
plt.xlim(left=-250, right=1250)
plt.xticks(ticks=ticks, labels=labels)
plt.ylabel('spike-in normalized dyad occupancy')
plt.ylim(0, 40)
plt.yticks(ticks=np.array((0, 10,20,30,40)))
plt.legend(loc='upper right')
plt.title('POB3')
plt.tight_layout()
plt.savefig('plots/metagenes/POB3_DMSO.png', dpi=300)
plt.show();

plt.figure(dpi=300)
plt.plot(positions, spt6_e154k_dmso_averaged.mean(), label = 'SPT6 DMSO', color = 'black', lw=2)
plt.plot(positions, dntd_e154k_dmso_averaged.mean(), label = 'ΔNTD DMSO', color = 'red', lw=2)
plt.plot(positions, ev_e154k_dmso_averaged.mean(), label = 'EV DMSO', color = 'darkorange', lw=2)
plt.xlabel('position')
plt.xlim(left=-250, right=1250)
plt.xticks(ticks=ticks, labels=labels)
plt.ylabel('spike-in normalized dyad occupancy')
plt.ylim(0, 40)
plt.yticks(ticks=np.array((0, 10,20,30,40)))
plt.legend(loc='upper right')
plt.title('E154K')
plt.tight_layout()
plt.savefig('plots/metagenes/E154K_DMSO.png', dpi=300)
plt.show();





plt.figure(dpi=300)
plt.plot(positions, spt6_pob3_iaa_averaged.mean(), label = 'SPT6 POB3', color = 'black', lw=2)
plt.plot(positions, spt6_e154k_iaa_averaged.mean(), label = 'SPT6 E154K', color = 'dodgerblue', lw=2)
plt.xlabel('position')
plt.xlim(left=-250, right=1250)
plt.xticks(ticks=ticks, labels=labels)
plt.ylabel('spike-in normalized dyad occupancy')
plt.ylim(0, 40)
plt.yticks(ticks=np.array((0, 10,20,30,40)))
plt.legend(loc='upper right')
plt.title('IAA')
plt.tight_layout()
plt.savefig('plots/metagenes/SPT6_IAA.png', dpi=300)
plt.show();

plt.figure(dpi=300)
plt.plot(positions, dntd_pob3_iaa_averaged.mean(), label = 'ΔNTD POB3', color = 'red', lw=2)
plt.plot(positions, dntd_e154k_iaa_averaged.mean(), label = 'ΔNTD E154K', color = 'orange', lw=2)
plt.xlabel('position')
plt.xlim(left=-250, right=1250)
plt.xticks(ticks=ticks, labels=labels)
plt.ylabel('spike-in normalized dyad occupancy')
plt.ylim(0, 40)
plt.yticks(ticks=np.array((0, 10,20,30,40)))
plt.legend(loc='upper right')
plt.title('IAA')
plt.tight_layout()
plt.savefig('plots/metagenes/DNTD_IAA.png', dpi=300)
plt.show();

plt.figure(dpi=300)
plt.plot(positions, ev_pob3_iaa_averaged.mean(), label = 'EV POB3', color = 'green', lw=2)
plt.plot(positions, ev_e154k_iaa_averaged.mean(), label = 'EV E154K', color = 'purple', lw=2)
plt.xlabel('position')
plt.xlim(left=-250, right=1250)
plt.xticks(ticks=ticks, labels=labels)
plt.ylabel('spike-in normalized dyad occupancy')
plt.ylim(0, 40)
plt.yticks(ticks=np.array((0, 10,20,30,40)))
plt.legend(loc='upper right')
plt.title('IAA')
plt.tight_layout()
plt.savefig('plots/metagenes/EV_IAA.png', dpi=300)
plt.show();





plt.figure(dpi=300)
plt.plot(positions, spt6_pob3_dmso_averaged.mean(), label = 'SPT6 POB3', color = 'black', lw=2)
plt.plot(positions, spt6_e154k_dmso_averaged.mean(), label = 'SPT6 E154K', color = 'dodgerblue', lw=2)
plt.xlabel('position')
plt.xlim(left=-250, right=1250)
plt.xticks(ticks=ticks, labels=labels)
plt.ylabel('spike-in normalized dyad occupancy')
plt.ylim(0, 40)
plt.yticks(ticks=np.array((0, 10,20,30,40)))
plt.legend(loc='upper right')
plt.title('DMSO')
plt.tight_layout()
plt.savefig('plots/metagenes/SPT6_DMSO.png', dpi=300)
plt.show();

plt.figure(dpi=300)
plt.plot(positions, dntd_pob3_dmso_averaged.mean(), label = 'ΔNTD POB3', color = 'red', lw=2)
plt.plot(positions, dntd_e154k_dmso_averaged.mean(), label = 'ΔNTD E154K', color = 'darkorange', lw=2)
plt.xlabel('position')
plt.xlim(left=-250, right=1250)
plt.xticks(ticks=ticks, labels=labels)
plt.ylabel('spike-in normalized dyad occupancy')
plt.ylim(0, 40)
plt.yticks(ticks=np.array((0, 10,20,30,40)))
plt.legend(loc='upper right')
plt.title('DMSO')
plt.tight_layout()
plt.savefig('plots/metagenes/DNTD_DMSO.png', dpi=300)
plt.show();

plt.figure(dpi=300)
plt.plot(positions, ev_pob3_dmso_averaged.mean(), label = 'EV POB3', color = 'green', lw=2)
plt.plot(positions, ev_e154k_dmso_averaged.mean(), label = 'EV E154K', color = 'purple', lw=2)
plt.xlabel('position')
plt.xlim(left=-250, right=1250)
plt.xticks(ticks=ticks, labels=labels)
plt.ylabel('spike-in normalized dyad occupancy')
plt.ylim(0, 40)
plt.yticks(ticks=np.array((0, 10,20,30,40)))
plt.legend(loc='upper right')
plt.title('DMSO')
plt.tight_layout()
plt.savefig('plots/metagenes/EV_DMSO.png', dpi=300)
plt.show();



#heatmaps
ticks = [0,50,100,151]
labels = ['-0.5 kb', '+1 dyad', '+0.5 kb', '+1 kb']
sns.set_style('ticks', {'axes.facecolor':'white'})

minimum = np.min([spt6_pob3_iaa_averaged.quantile([0.10]).min(1),
                  ev_pob3_iaa_averaged.quantile([0.10]).min(1),
                  dntd_pob3_iaa_averaged.quantile([0.10]).min(1),
                  spt6_e154k_iaa_averaged.quantile([0.10]).min(1),
                  ev_e154k_iaa_averaged.quantile([0.10]).min(1),
                  dntd_e154k_iaa_averaged.quantile([0.10]).min(1)
                  ])
maximum = np.max([spt6_pob3_iaa_averaged.quantile([0.85]).max(1),
                  ev_pob3_iaa_averaged.quantile([0.85]).max(1),
                  dntd_pob3_iaa_averaged.quantile([0.85]).max(1),
                  spt6_e154k_iaa_averaged.quantile([0.85]).max(1),
                  ev_e154k_iaa_averaged.quantile([0.85]).max(1),
                  dntd_e154k_iaa_averaged.quantile([0.85]).max(1)
                  ])


minimum = 0
maximum = 60

#colorbar
a = np.array([[0,60]])
plt.figure(figsize=(1,4), dpi=600)
img = plt.imshow(a, cmap="afmhot")
plt.gca().set_visible(False)
cbar = plt.colorbar(ticks=[0,20,40,60],
                    orientation='vertical')
cbar.ax.set_yticklabels(['$\\leq$$0$', '20', '40', '$\\geq$$60$'], fontsize=6)
plt.tight_layout()
plt.savefig('plots/heatmaps/heatmaps_colorbar.png', dpi=600)
plt.show();

#SPT6 POB3
f = plt.figure(figsize=(2.8,4), dpi=600)
gs = f.add_gridspec(1,1)

ax = f.add_subplot(gs[0,0])
sns.heatmap(data=spt6_pob3_iaa_averaged.iloc[:,0:151], cmap=sns.color_palette("afmhot", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=0.5, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=10)

f.tight_layout()
plt.savefig("plots/heatmaps/spt6_pob3_Rpb1sorted.png", dpi=600)
plt.show();


#Δ2-238 POB3
f = plt.figure(figsize=(2.8,4), dpi=600)
gs = f.add_gridspec(1,1)

ax = f.add_subplot(gs[0,0])
sns.heatmap(data=dntd_pob3_iaa_averaged.iloc[:,0:151], cmap=sns.color_palette("afmhot", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=0.5, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=10)

f.tight_layout()
plt.savefig("plots/heatmaps/dntd_pob3_Rpb1sorted.png", dpi=600)
plt.show();

#EV POB3
f = plt.figure(figsize=(2.8,4), dpi=600)
gs = f.add_gridspec(1,1)

ax = f.add_subplot(gs[0,0])
sns.heatmap(data=ev_pob3_iaa_averaged.iloc[:,0:151], cmap=sns.color_palette("afmhot", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=0.5, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=10)

f.tight_layout()
plt.savefig("plots/heatmaps/ev_pob3_Rpb1sorted.png", dpi=600)
plt.show();

#SPT6 E154K
f = plt.figure(figsize=(2.8,4), dpi=600)
gs = f.add_gridspec(1,1)

ax = f.add_subplot(gs[0,0])
sns.heatmap(data=spt6_e154k_iaa_averaged.iloc[:,0:151], cmap=sns.color_palette("afmhot", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=0.5, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=10)

f.tight_layout()
plt.savefig("plots/heatmaps/spt6_e154k_Rpb1sorted.png", dpi=600)
plt.show();

#Δ2-238 E154K
f = plt.figure(figsize=(2.8,4), dpi=600)
gs = f.add_gridspec(1,1)

ax = f.add_subplot(gs[0,0])
sns.heatmap(data=dntd_e154k_iaa_averaged.iloc[:,0:151], cmap=sns.color_palette("afmhot", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=0.5, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=10)

f.tight_layout()
plt.savefig("plots/heatmaps/dntd_e154k_Rpb1sorted.png", dpi=600)
plt.show();

#EV E154K
f = plt.figure(figsize=(2.8,4), dpi=600)
gs = f.add_gridspec(1,1)

ax = f.add_subplot(gs[0,0])
sns.heatmap(data=ev_e154k_iaa_averaged.iloc[:,0:151], cmap=sns.color_palette("afmhot", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=0.5, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=10)

f.tight_layout()
plt.savefig("plots/heatmaps/ev_e154k_Rpb1sorted.png", dpi=600)
plt.show();






##Other heatmaps
dmso_minimum = np.min([spt6_pob3_dmso_averaged.quantile([0.10]).min(1),
                  ev_pob3_dmso_averaged.quantile([0.10]).min(1),
                  dntd_pob3_dmso_averaged.quantile([0.10]).min(1)
                  ])
dmso_maximum = np.max([spt6_pob3_dmso_averaged.quantile([0.85]).max(1),
                  ev_pob3_dmso_averaged.quantile([0.85]).max(1),
                  dntd_pob3_dmso_averaged.quantile([0.85]).max(1)
                  ])

iaa_minimum = np.min([spt6_pob3_iaa_averaged.quantile([0.10]).min(1),
                  ev_pob3_iaa_averaged.quantile([0.10]).min(1),
                  dntd_pob3_iaa_averaged.quantile([0.10]).min(1)
                  ])
iaa_maximum = np.max([spt6_pob3_iaa_averaged.quantile([0.85]).max(1),
                  ev_pob3_iaa_averaged.quantile([0.85]).max(1),
                  dntd_pob3_iaa_averaged.quantile([0.85]).max(1)
                  ])

minimum = np.min([dmso_minimum,iaa_minimum])
maximum = np.max([dmso_maximum,iaa_maximum])

f = plt.figure(figsize=(8.5,11), dpi=600)
gs = f.add_gridspec(2,3)
f.suptitle('IAA', fontsize = 12)

ax = f.add_subplot(gs[0,0])
sns.heatmap(data=spt6_pob3_iaa_averaged.iloc[:,0:151], cmap=sns.color_palette("afmhot", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=0.5, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=10)
plt.title('SPT6 POB3', fontsize=10)    

ax = f.add_subplot(gs[0,1])
sns.heatmap(data=dntd_pob3_iaa_averaged.iloc[:,0:151], cmap=sns.color_palette("afmhot", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=0.5, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=10)
plt.title('ΔNTD POB3', fontsize=10)

ax = f.add_subplot(gs[0,2])
sns.heatmap(data=ev_pob3_iaa_averaged.iloc[:,0:151], cmap=sns.color_palette("afmhot", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=0.5, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=10)
plt.title('EV POB3', fontsize=10)    

ax = f.add_subplot(gs[1,0])
sns.heatmap(data=spt6_e154k_iaa_averaged.iloc[:,0:151], cmap=sns.color_palette("afmhot", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=0.5, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=10)
plt.title('SPT6 E154K', fontsize=10)    

ax = f.add_subplot(gs[1,1])
sns.heatmap(data=dntd_e154k_iaa_averaged.iloc[:,0:151], cmap=sns.color_palette("afmhot", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=0.5, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=10)
plt.title('ΔNTD E154K', fontsize=10)

ax = f.add_subplot(gs[1,2])
sns.heatmap(data=ev_e154k_iaa_averaged.iloc[:,0:151], cmap=sns.color_palette("afmhot", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=0.5, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=10)
plt.title('EV E154K', fontsize=10)    

f.tight_layout()
plt.savefig("plots/heatmaps/IAA_Rpb1sorted.png", dpi=600)
plt.show();




f = plt.figure(figsize=(8.5,11), dpi=600)
gs = f.add_gridspec(2,3)
f.suptitle('DMSO', fontsize = 12)

ax = f.add_subplot(gs[0,0])
sns.heatmap(data=spt6_pob3_dmso_averaged.iloc[:,0:151], cmap=sns.color_palette("afmhot", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=0.5, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=10)
plt.title('SPT6 POB3', fontsize=10)    

ax = f.add_subplot(gs[0,1])
sns.heatmap(data=dntd_pob3_dmso_averaged.iloc[:,0:151], cmap=sns.color_palette("afmhot", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=0.5, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=10)
plt.title('ΔNTD POB3', fontsize=10)

ax = f.add_subplot(gs[0,2])
sns.heatmap(data=ev_pob3_dmso_averaged.iloc[:,0:151], cmap=sns.color_palette("afmhot", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=0.5, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=10)
plt.title('EV POB3', fontsize=10)    

ax = f.add_subplot(gs[1,0])
sns.heatmap(data=spt6_e154k_dmso_averaged.iloc[:,0:151], cmap=sns.color_palette("afmhot", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=0.5, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=10)
plt.title('SPT6 E154K', fontsize=10)    

ax = f.add_subplot(gs[1,1])
sns.heatmap(data=dntd_e154k_dmso_averaged.iloc[:,0:151], cmap=sns.color_palette("afmhot", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=0.5, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=10)
plt.title('ΔNTD E154K', fontsize=10)

ax = f.add_subplot(gs[1,2])
sns.heatmap(data=ev_e154k_dmso_averaged.iloc[:,0:151], cmap=sns.color_palette("afmhot", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=0.5, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=10)
plt.title('EV E154K', fontsize=10)    

f.tight_layout()
plt.savefig("plots/heatmaps/DMSO_Rpb1sorted.png", dpi=600)
plt.show();




#fold-change
minimum = -2
maximum = 2

f = plt.figure(figsize=(8.5,11), dpi=600)
gs = f.add_gridspec(2,3)
f.suptitle('IAA', fontsize = 12)

ax = f.add_subplot(gs[0,0])
sns.heatmap(data=np.log2((dntd_pob3_iaa_averaged.iloc[:,0:151]+0.01)/(spt6_pob3_iaa_averaged.iloc[:,0:151]+0.01)),
            cmap=sns.color_palette("bwr", as_cmap=True),
            vmax=maximum,
            vmin=minimum,
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=0.5, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=8)
plt.title('ΔNTD/WT', fontsize=10)

ax = f.add_subplot(gs[0,1])
sns.heatmap(data=np.log2((ev_pob3_iaa_averaged.iloc[:,0:151]+0.01)/(spt6_pob3_iaa_averaged.iloc[:,0:151]+0.01)),
            cmap=sns.color_palette("bwr", as_cmap=True),
            vmax=maximum,
            vmin=minimum,
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=0.5, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=8)
plt.title('EV/WT', fontsize=10)

ax = f.add_subplot(gs[1,0])
sns.heatmap(data=np.log2((dntd_e154k_iaa_averaged.iloc[:,0:151]+0.01)/(spt6_e154k_iaa_averaged.iloc[:,0:151]+0.01)),
            cmap=sns.color_palette("bwr", as_cmap=True),
            vmax=maximum,
            vmin=minimum,
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=0.5, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=8)
plt.title('ΔNTD E154K/E154K', fontsize=10)


ax = f.add_subplot(gs[1,1])
sns.heatmap(data=np.log2((ev_e154k_iaa_averaged.iloc[:,0:151]+0.01)/(spt6_e154k_iaa_averaged.iloc[:,0:151]+0.01)),
            cmap=sns.color_palette("bwr", as_cmap=True),
            vmax=maximum,
            vmin=minimum,
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=0.5, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=8)
plt.title('EV E154K/E154K', fontsize=10)

ax = f.add_subplot(gs[0,2])
sns.heatmap(data=np.log2((spt6_e154k_iaa_averaged.iloc[:,0:151]+0.01)/(spt6_pob3_iaa_averaged.iloc[:,0:151]+0.01)),
            cmap=sns.color_palette("bwr", as_cmap=True),
            vmax=maximum,
            vmin=minimum,
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=0.5, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=8)
plt.title('E154K/WT', fontsize=10)


f.tight_layout()
plt.savefig("plots/heatmaps/logfc_IAA.png", dpi=600)
plt.show();


fc = np.log2((spt6_e154k_iaa_averaged.iloc[:,0:151]+0.01)/(spt6_pob3_iaa_averaged.iloc[:,0:151]+0.01))
f = plt.figure(figsize=(6,8), dpi=600)
gs = f.add_gridspec(1,1)

ax = f.add_subplot(gs[0,0])
sns.heatmap(data=fc, cmap=sns.color_palette("bwr", as_cmap=True),
            vmax=maximum,
            vmin=minimum,
            yticklabels=False, cbar=True)

for spine in ax.spines.values():
    spine.set(visible=True, lw=0.5, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=8)
plt.title('E154K/WT', fontsize=22)

f.tight_layout()
#plt.savefig("plots/heatmaps/reference/lala.png", dpi=600)
plt.show();





#Sorted by gene length
positions = np.arange(-500,4500,10)

spt6_pob3_iaa_averaged = pd.read_csv('/Users/jlwarner/Desktop/mnase_seq_new/transfer/tab/425_93_I_averaged_plusone_reference.tab', 
            sep='\t', header=3, names=positions)
ev_pob3_iaa_averaged = pd.read_csv('/Users/jlwarner/Desktop/mnase_seq_new/transfer/tab/425_414_I_averaged_plusone_reference.tab', 
            sep='\t', header=3, names=positions)
dntd_pob3_iaa_averaged = pd.read_csv('/Users/jlwarner/Desktop/mnase_seq_new/transfer/tab/425_95_I_averaged_plusone_reference.tab', 
            sep='\t', header=3, names=positions)
spt6_e154k_iaa_averaged = pd.read_csv('/Users/jlwarner/Desktop/mnase_seq_new/transfer/tab/553_93_I_averaged_plusone_reference.tab', 
            sep='\t', header=3, names=positions)
ev_e154k_iaa_averaged = pd.read_csv('/Users/jlwarner/Desktop/mnase_seq_new/transfer/tab/553_414_I_averaged_plusone_reference.tab', 
            sep='\t', header=3, names=positions)
dntd_e154k_iaa_averaged = pd.read_csv('/Users/jlwarner/Desktop/mnase_seq_new/transfer/tab/553_95_I_averaged_plusone_reference.tab', 
            sep='\t', header=3, names=positions)

spt6_pob3_dmso_averaged = pd.read_csv('/Users/jlwarner/Desktop/mnase_seq_new/transfer/tab/425_93_D_averaged_plusone_reference.tab', 
            sep='\t', header=3, names=positions)
ev_pob3_dmso_averaged = pd.read_csv('/Users/jlwarner/Desktop/mnase_seq_new/transfer/tab/425_414_D_averaged_plusone_reference.tab', 
            sep='\t', header=3, names=positions)
dntd_pob3_dmso_averaged = pd.read_csv('/Users/jlwarner/Desktop/mnase_seq_new/transfer/tab/425_95_D_averaged_plusone_reference.tab', 
            sep='\t', header=3, names=positions)
spt6_e154k_dmso_averaged = pd.read_csv('/Users/jlwarner/Desktop/mnase_seq_new/transfer/tab/553_93_D_averaged_plusone_reference.tab', 
            sep='\t', header=3, names=positions)
ev_e154k_dmso_averaged = pd.read_csv('/Users/jlwarner/Desktop/mnase_seq_new/transfer/tab/553_414_D_averaged_plusone_reference.tab', 
            sep='\t', header=3, names=positions)
dntd_e154k_dmso_averaged = pd.read_csv('/Users/jlwarner/Desktop/mnase_seq_new/transfer/tab/553_95_D_averaged_plusone_reference.tab', 
            sep='\t', header=3, names=positions)

ticks = [0,50,150,250,350,450]
labels = ['-0.5 kb', '+1 dyad', '+1 kb','+2 kb', '+3 kb', '+4 kb']
sns.set_style('ticks', {'axes.facecolor':'#DCDCDC'})

dmso_minimum = np.min([spt6_pob3_dmso_averaged.quantile([0.05]).min(1),
                  ev_pob3_dmso_averaged.quantile([0.05]).min(1),
                  dntd_pob3_dmso_averaged.quantile([0.05]).min(1)
                  ])
dmso_maximum = np.max([spt6_pob3_dmso_averaged.quantile([0.85]).max(1),
                  ev_pob3_dmso_averaged.quantile([0.85]).max(1),
                  dntd_pob3_dmso_averaged.quantile([0.85]).max(1)
                  ])

iaa_minimum = np.min([spt6_pob3_iaa_averaged.quantile([0.05]).min(1),
                  ev_pob3_iaa_averaged.quantile([0.05]).min(1),
                  dntd_pob3_iaa_averaged.quantile([0.05]).min(1)
                  ])
iaa_maximum = np.max([spt6_pob3_iaa_averaged.quantile([0.85]).max(1),
                  ev_pob3_iaa_averaged.quantile([0.85]).max(1),
                  dntd_pob3_iaa_averaged.quantile([0.85]).max(1)
                  ])

minimum = np.min([dmso_minimum,iaa_minimum])
maximum = np.max([dmso_maximum,iaa_maximum])

f = plt.figure(figsize=(8.5,11), dpi=600)
gs = f.add_gridspec(2,3)
f.suptitle('IAA', fontsize = 12)

ax = f.add_subplot(gs[0,0])
sns.heatmap(data=spt6_pob3_iaa_averaged, cmap=sns.color_palette("afmhot", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=0.5, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=10)
plt.title('SPT6 POB3', fontsize=10)    

ax = f.add_subplot(gs[0,1])
sns.heatmap(data=dntd_pob3_iaa_averaged, cmap=sns.color_palette("afmhot", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=0.5, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=10)
plt.title('ΔNTD POB3', fontsize=10)

ax = f.add_subplot(gs[0,2])
sns.heatmap(data=ev_pob3_iaa_averaged, cmap=sns.color_palette("afmhot", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=0.5, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=10)
plt.title('EV POB3', fontsize=10)    

ax = f.add_subplot(gs[1,0])
sns.heatmap(data=spt6_e154k_iaa_averaged, cmap=sns.color_palette("afmhot", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=0.5, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=10)
plt.title('SPT6 E154K', fontsize=10)    

ax = f.add_subplot(gs[1,1])
sns.heatmap(data=dntd_e154k_iaa_averaged, cmap=sns.color_palette("afmhot", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=0.5, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=10)
plt.title('ΔNTD E154K', fontsize=10)

ax = f.add_subplot(gs[1,2])
sns.heatmap(data=ev_e154k_iaa_averaged, cmap=sns.color_palette("afmhot", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=0.5, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=10)
plt.title('EV E154K', fontsize=10)    

f.tight_layout()
#plt.savefig("plots/heatmaps/IAA_lengthsorted.png", dpi=600)
plt.show();




f = plt.figure(figsize=(8.5,11), dpi=600)
gs = f.add_gridspec(2,3)
f.suptitle('DMSO', fontsize = 12)

ax = f.add_subplot(gs[0,0])
sns.heatmap(data=spt6_pob3_dmso_averaged, cmap=sns.color_palette("afmhot", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=0.5, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=10)
plt.title('SPT6 POB3', fontsize=10)    

ax = f.add_subplot(gs[0,1])
sns.heatmap(data=dntd_pob3_dmso_averaged, cmap=sns.color_palette("afmhot", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=0.5, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=10)
plt.title('ΔNTD POB3', fontsize=10)

ax = f.add_subplot(gs[0,2])
sns.heatmap(data=ev_pob3_dmso_averaged, cmap=sns.color_palette("afmhot", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=0.5, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=10)
plt.title('EV POB3', fontsize=10)    

ax = f.add_subplot(gs[1,0])
sns.heatmap(data=spt6_e154k_dmso_averaged, cmap=sns.color_palette("afmhot", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=0.5, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=10)
plt.title('SPT6 E154K', fontsize=10)    

ax = f.add_subplot(gs[1,1])
sns.heatmap(data=dntd_e154k_dmso_averaged, cmap=sns.color_palette("afmhot", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=0.5, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=10)
plt.title('ΔNTD E154K', fontsize=10)

ax = f.add_subplot(gs[1,2])
sns.heatmap(data=ev_e154k_dmso_averaged, cmap=sns.color_palette("afmhot", as_cmap=True), 
            vmax=maximum,
            vmin=minimum, 
            yticklabels=False, cbar=False)

for spine in ax.spines.values():
    spine.set(visible=True, lw=0.5, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=10)
plt.title('EV E154K', fontsize=10)    

f.tight_layout()
#plt.savefig("plots/heatmaps/DMSO_lengthsorted.png", dpi=600)
plt.show();



#fold-change
fc = np.log2(ev_pob3_iaa_averaged/spt6_pob3_iaa_averaged)
minimum = -3
maximum = 3

f = plt.figure(figsize=(6,8), dpi=600)
gs = f.add_gridspec(1,1)
#f.suptitle('log2 $\\frac{ΔNTD\:IAA}{SPT6\:IAA}$ Rpb1 ChIP occupancy', fontsize = 30)

ax = f.add_subplot(gs[0,0])
sns.heatmap(data=fc, cmap=sns.color_palette("bwr", as_cmap=True),
            vmax=maximum,
            vmin=minimum,
            yticklabels=False, cbar=True)

for spine in ax.spines.values():
    spine.set(visible=True, lw=0.5, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=8)
plt.title('EV/WT', fontsize=22)

f.tight_layout()
#plt.savefig("plots/heatmaps/reference/lala.png", dpi=600)
plt.show();

#fold-change
fc = np.log2(dntd_pob3_iaa_averaged/spt6_pob3_iaa_averaged)
minimum = -3
maximum = 3

f = plt.figure(figsize=(6,8), dpi=600)
gs = f.add_gridspec(1,1)
#f.suptitle('log2 $\\frac{ΔNTD\:IAA}{SPT6\:IAA}$ Rpb1 ChIP occupancy', fontsize = 30)

ax = f.add_subplot(gs[0,0])
sns.heatmap(data=fc, cmap=sns.color_palette("bwr", as_cmap=True),
            vmax=maximum,
            vmin=minimum,
            yticklabels=False, cbar=True)

for spine in ax.spines.values():
    spine.set(visible=True, lw=0.5, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=8)
plt.title('ΔNTD/WT', fontsize=22)

f.tight_layout()
#plt.savefig("plots/heatmaps/reference/lala.png", dpi=600)
plt.show();

#fold-change
fc = np.log2(dntd_e154k_iaa_averaged/dntd_pob3_iaa_averaged)
minimum = -3
maximum = 3

f = plt.figure(figsize=(6,8), dpi=600)
gs = f.add_gridspec(1,1)
#f.suptitle('log2 $\\frac{ΔNTD\:IAA}{SPT6\:IAA}$ Rpb1 ChIP occupancy', fontsize = 30)

ax = f.add_subplot(gs[0,0])
sns.heatmap(data=fc, cmap=sns.color_palette("bwr", as_cmap=True),
            vmax=maximum,
            vmin=minimum,
            yticklabels=False, cbar=True)

for spine in ax.spines.values():
    spine.set(visible=True, lw=0.5, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=8)
plt.title('ΔNTD E154K/ΔNTD', fontsize=22)

f.tight_layout()
#plt.savefig("plots/heatmaps/reference/lala.png", dpi=600)
plt.show();

#fold-change
fc = np.log2(dntd_e154k_iaa_averaged/spt6_pob3_iaa_averaged)
minimum = -3
maximum = 3

f = plt.figure(figsize=(6,8), dpi=600)
gs = f.add_gridspec(1,1)
#f.suptitle('log2 $\\frac{ΔNTD\:IAA}{SPT6\:IAA}$ Rpb1 ChIP occupancy', fontsize = 30)

ax = f.add_subplot(gs[0,0])
sns.heatmap(data=fc, cmap=sns.color_palette("bwr", as_cmap=True),
            vmax=maximum,
            vmin=minimum,
            yticklabels=False, cbar=True)

for spine in ax.spines.values():
    spine.set(visible=True, lw=0.5, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=8)
plt.title('ΔNTD E154K/WT', fontsize=22)

f.tight_layout()
#plt.savefig("plots/heatmaps/reference/lala.png", dpi=600)
plt.show();

#fold-change
fc = np.log2(spt6_e154k_iaa_averaged/spt6_pob3_iaa_averaged)
minimum = -3
maximum = 3

f = plt.figure(figsize=(6,8), dpi=600)
gs = f.add_gridspec(1,1)
#f.suptitle('log2 $\\frac{ΔNTD\:IAA}{SPT6\:IAA}$ Rpb1 ChIP occupancy', fontsize = 30)

ax = f.add_subplot(gs[0,0])
sns.heatmap(data=fc, cmap=sns.color_palette("bwr", as_cmap=True),
            vmax=maximum,
            vmin=minimum,
            yticklabels=False, cbar=True)

for spine in ax.spines.values():
    spine.set(visible=True, lw=0.5, color='black')
#plt.ylabel('genes', fontsize=20)
plt.xticks(ticks=ticks, labels=labels, fontsize=8)
plt.title('E154K/WT', fontsize=22)

f.tight_layout()
#plt.savefig("plots/heatmaps/reference/lala.png", dpi=600)
plt.show();



