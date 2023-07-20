#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# Calculates Average normalized H3K4me3 for ILC2s generated for this manuscript at 
# TSSs specific to ILC1s and ILC2s.
#
# Generates Figure 1c plot.


# In[ ]:


import pandas as pd
import glob
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
import seaborn as sns

# deeptools 3.5.1 used for this analysis 


# ### Output Directory

# In[ ]:


# Create Output Dir variable
output_dir = './Figure_1/'


# ### Compare in vitro ILC2 H3K4me3 at ILC1 and ILC2 specific marked TSSs
# 1. in vitro ILC2 data were generated for this manuscript

# In[ ]:


# Store paths to ILC2 H3K4me3 bigWig files
K4me3_bw = sorted(glob.glob('./processed_data/blfH3K4ME3_rep*bw'))

# Convert bigWig file name list to space delimited string
str_K4me3_bw = " ".join(K4me3_bw)

# Experiment names 
K4me3_exp_names = ["H3K4ME3_rep1", "H3K4ME3_rep2"]

# Convert experiment names list to space delimited string
str_K4me3_exp_names = " ".join(K4me3_exp_names)


# In[ ]:


# Differential TSSs marked H3K4me3 specific to ILC1s 
ILC1_TSS_BED = f'{output_dir}ILC1_specific_k4me3-marked_TSSs.bed'


# In[ ]:


get_ipython().run_cell_magic('bash', '-s  "$str_K4me3_bw" "$str_K4me3_exp_names" "$ILC1_TSS_BED" "$output_dir"', '\n# Calculate normalized H3K4me3 signal at TSSs specific to ILC1s\nmultiBigwigSummary BED-file \\\n--bwfiles $1 \\\n--BED $3 \\\n--labels $2 \\\n--outFileName $4ILC2_H3K4me3_at_GBA_K4me3_marked_TSSs_specificToILC1.npz \\\n--outRawCounts $4ILC2_H3K4me3_at_GBA_K4me3_marked_TSSs_specificToILC1.tab \\\n-p 4\n')


# In[ ]:


# Differential TSSs marked H3K4me3 specific to ILC2s 
ILC2_TSS_BED = f'{output_dir}ILC2_specific_k4me3-marked_TSSs.bed'


# In[ ]:


get_ipython().run_cell_magic('bash', '-s  "$str_K4me3_bw" "$str_K4me3_exp_names" "$ILC2_TSS_BED" "$output_dir"', '\n# Calculate normalized H3K4me3 signal at TSSs specific to ILC2s\nmultiBigwigSummary BED-file \\\n--bwfiles $1 \\\n--BED $3 \\\n--labels $2 \\\n--outFileName $4ILC2_H3K4me3_at_GBA_K4me3_marked_TSSs_specificToILC2.npz \\\n--outRawCounts $4ILC2_H3K4me3_at_GBA_K4me3_marked_TSSs_specificToILC2.tab \\\n-p 4\n')


# In[ ]:


# Read ILC2 normalized H3K4me3 signal at ILC1 TSSs 
K4me3_ILC1 = pd.read_csv(f'{output_dir}ILC2_H3K4me3_at_GBA_K4me3_marked_TSSs_specificToILC1.tab',
                                sep = '\t',
                               header = 0,
                               names = ['chr','start','end','H3K4me3_ILC1_rep1','H3K4me3_ILC1_rep2'])

# Average normalized H3K4me3 signal across both replicates
K4me3_ILC1 = K4me3_ILC1.assign(H3K4me3_ILC1_avg = K4me3_ILC1.loc[:,'H3K4me3_ILC1_rep1':"H3K4me3_ILC1_rep2"].mean(axis = 1))

# Drop the individual replicate data
K4me3_ILC1 = K4me3_ILC1.drop(columns=['H3K4me3_ILC1_rep1',"H3K4me3_ILC1_rep2"])

# Wide to Long DF
K4me3_ILC1_long = pd.melt(K4me3_ILC1, id_vars=['chr','start','end'], var_name='cellType', value_name= 'k4me3_signal')


# In[ ]:


# Read ILC2 normalized H3K4me3 signal at ILC2 TSSs 
K4me3_ILC2 = pd.read_csv(f'{output_dir}ILC2_H3K4me3_at_GBA_K4me3_marked_TSSs_specificToILC2.tab',
                                sep = '\t',
                               header = 0,
                               names = ['chr','start','end','H3K4me3_ILC2_rep1','H3K4me3_ILC2_rep2'])

# Average normalized H3K4me3 signal across both replicates
K4me3_ILC2 = K4me3_ILC2.assign(H3K4me3_ILC2_avg = K4me3_ILC2.loc[:,'H3K4me3_ILC2_rep1':"H3K4me3_ILC2_rep2"].mean(axis = 1))

# Drop the individual replicate data
K4me3_ILC2 = K4me3_ILC2.drop(columns=['H3K4me3_ILC2_rep1',"H3K4me3_ILC2_rep2"])

# Wide to Long DF
K4me3_ILC2_long = pd.melt(K4me3_ILC2, id_vars=['chr','start','end'], var_name='cellType', value_name= 'k4me3_signal')


# In[ ]:


# Concat DFs
K4me3_ILC1_ILC2_long = pd.concat([K4me3_ILC1_long,K4me3_ILC2_long])


# In[ ]:


fix, ax = plt.subplots(figsize=(2.5,3))


sns.boxplot(x="cellType",
            y="k4me3_signal",
            data=K4me3_ILC1_ILC2_long,
            ax = ax,
           showfliers = False,
           linewidth = 1,
           width = .5)

ax.set_xlabel("")
ax.set_ylabel("H3K4me3 \n (Normalized signal)")
ax.set_xticklabels(["ILC1", "ILC2"])

sns.despine()

plt.savefig(f'{output_dir}Fig1c_ILC2_H3K4me3_at_GBA_ILC_specific_TSSs.pdf',
            bbox_inches="tight",
            transparent=True)


# In[ ]:


#Calculate Mann-Whitney rank test 
from scipy.stats import mannwhitneyu

mannwhitneyu(K4me3_ILC1['H3K4me3_ILC1_avg'],
             K4me3_ILC2['H3K4me3_ILC2_avg'],
            alternative = "less")



# In[ ]:




