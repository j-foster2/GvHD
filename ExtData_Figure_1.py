#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# Performs hierarchical clustering on H3K4me3 signal at differential TSS across ILCs
#
# Generates plot for Extended Data Figure 1c-e


# In[ ]:


import pandas as pd
import numpy as np
import glob
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
import seaborn as sns
from scipy.cluster import hierarchy
from scipy.cluster.hierarchy import dendrogram
from sklearn import preprocessing
from scipy.cluster.hierarchy import fcluster
from distinctipy import distinctipy
from scipy.stats import mannwhitneyu


# ### Output Directory

# In[ ]:


# Create Output Dir variable
output_dir = './Figure_1/'


# ### Generate BED file of TSSs with differential H3K4me3 across ILCs

# In[ ]:


# Read ILC1 DESeq2 output
ILC_LRT_deseq2 = pd.read_csv(f'{output_dir}ILC_LRT_deseq2_results.txt',
                         sep = '\t')

# Format gene name column
ILC_LRT_deseq2 = ILC_LRT_deseq2.rename(columns = {"Unnamed: 0":"geneName"})

# Set index to gene names
ILC_LRT_deseq2 = ILC_LRT_deseq2.set_index("geneName")

# Filter for TSSs with significantly differential H3K4me3 signal across ILCs (GBA) 
ILC_LRT_deseq2 = ILC_LRT_deseq2[(ILC_LRT_deseq2['padj'] <= 0.05)].sort_values(by = 'baseMean', ascending = False)

# Read RefSeq BED of all TSSs (including 300 bp upstream and 500 downstream)
refseq_bed = pd.read_csv('./Figure_1/refseq_promotor.bed',\
                        sep='\t',
                        header=0,
                        names = ['chr','start','end','geneName'],
                        index_col = 'geneName')

# Merge RefSeq BED and DESeq2 results
ILC_LRT_deseq2 = ILC_LRT_deseq2.join(refseq_bed)

# Write BED file of differnetial TSSs
ILC_LRT_deseq2.reset_index().to_csv(f'{output_dir}ILC_LRT_significant_TSSs.bed',
                   sep = '\t',
                   columns = ['chr','start','end','geneName'],
                  index = False,
                  header = None)


# ### Calculate H3K4me3 Signal at TSSs with differential H3K4me3 across ILCs

# In[ ]:


# Store bigWig file names for all GBA ILC H3K4me3 ChIP-seq data
gba_k4me3_bw = sorted(glob.glob('./processed_data/blfILC*GBA*bw'))

#Convert bigWig file name list to space delimited string
str_gba_k4me3_bw = " ".join(gba_k4me3_bw)

#Gather experiment names 
exp_names = ["ILC1_GBA_H3K4me3_rep1", "ILC1_GBA_H3K4me3_rep2",
             "ILC2_GBA_H3K4me3_rep1", "ILC2_GBA_H3K4me3_rep2",
             "ILC3_GBA_H3K4me3_rep1", "ILC3_GBA_H3K4me3_rep2"]

#Convert experiment names list to space delimited string
str_exp_names = " ".join(exp_names)


# In[ ]:


# Store path to BED of TSSs with differential H3K4me3 across ILCs 
ILC_LRT_TSS_bed = f'{output_dir}ILC_LRT_significant_TSSs.bed'


# In[ ]:


get_ipython().run_cell_magic('bash', '-s  "$str_gba_k4me3_bw" "$str_exp_names" "$ILC_LRT_TSS_bed" "$output_dir"', '\nmultiBigwigSummary BED-file \\\n--bwfiles $1 \\\n--BED $3 \\\n--labels $2 \\\n--outFileName $4GBA_ILC_H3K4me3_sig_TSSs.npz \\\n--outRawCounts $4GBA_ILC2_H3K4me3_sig_TSSs.tab \\\n-p 4\n')


# ### Identification of ILC1 and ILC2 specific k4me3-marked TSS

# In[ ]:


ILC_K4me3_TSS_data = pd.read_csv(f'{output_dir}GBA_ILC2_H3K4me3_sig_TSSs.tab', 
            sep = '\t',
            names = exp_names,
           header = 0).reset_index(drop = True)


# In[ ]:


# Calculate linkage matrix
K4me3_scaled_data = preprocessing.scale(ILC_K4me3_TSS_data)

Z = hierarchy.linkage(K4me3_scaled_data, method='average')

# Plot 
plt.figure()

plt.title("Dendrogram")

# Dendrogram plotting using linkage matrix
dendrogram2 = hierarchy.dendrogram(Z)


# In[ ]:


# Source of fancy_dendrogram function 
# https://joernhees.de/blog/2015/08/26/scipy-hierarchical-clustering-and-dendrogram-tutorial/

def fancy_dendrogram(*args, **kwargs):
    max_d = kwargs.pop('max_d', None)
    if max_d and 'color_threshold' not in kwargs:
        kwargs['color_threshold'] = max_d
    annotate_above = kwargs.pop('annotate_above', 0)

    ddata = dendrogram(*args, **kwargs)

    if not kwargs.get('no_plot', False):
        plt.title('Hierarchical Clustering Dendrogram (truncated)')
        plt.xlabel('sample index or (cluster size)')
        plt.ylabel('distance')
        for i, d, c in zip(ddata['icoord'], ddata['dcoord'], ddata['color_list']):
            x = 0.5 * sum(i[1:3])
            y = d[1]
            if y > annotate_above:
                plt.plot(x, y, 'o', c=c)
                plt.annotate("%.3g" % y, (x, y), xytext=(0, -5),
                             textcoords='offset points',
                             va='top', ha='center')
        if max_d:
            plt.axhline(y=max_d, c='k')
    return ddata


# In[ ]:


# Set Parameter to cut dendrogram
max_d = 1.7


# In[ ]:


fancy_dendrogram(
    Z,
    truncate_mode='lastp',
    p=12,
    leaf_rotation=90.,
    leaf_font_size=12.,
    show_contracted=True,
    annotate_above=10,
    max_d=max_d,  # plot a horizontal cut-off line
)
plt.show()


# In[ ]:


# Cluster asssignment for each differential TSS
clusters = fcluster(Z, max_d, criterion='distance')

# Read normalized H3K4me3 at differential TSSs
ILC_K4me3_TSS_data_clust = pd.read_csv(f'{output_dir}GBA_ILC2_H3K4me3_sig_TSSs.tab', 
            sep = '\t',
            names = exp_names,
           header = 0)

# Assign each differential TSSs to a cluster
ILC_K4me3_TSS_data_clust = ILC_K4me3_TSS_data_clust.assign(clusters = clusters)


# In[ ]:


# Unique Cluster names
cluster_names = np.unique(clusters)

#Isolate the promoters that are assigned to each cluster
cluster_df = [ILC_K4me3_TSS_data_clust[ILC_K4me3_TSS_data_clust['clusters'] == i] for i in cluster_names]

# Wide to Long DF conversion 
cluster_df_long = [pd.melt(df, id_vars=['clusters'], var_name='cellType', value_name= 'k4me3_signal') for df in cluster_df]


# #### Heatmap of H3K4me3 Clusters

# In[ ]:


cluster_data = pd.concat(cluster_df)

# Isolate clusters
clusters = cluster_data.pop("clusters")

# number of colours to generate
N = len(clusters.unique())

# generate N visually distinct colours
colors = distinctipy.get_colors(N, pastel_factor=1)


# In[ ]:


# Map Colors to clusters in Dataframe
lut = dict(zip(clusters.unique(), colors))

row_colors = clusters.map(lut)


# In[ ]:


g = sns.clustermap(cluster_data.loc[:,"ILC1_GBA_H3K4me3_rep1":"ILC2_GBA_H3K4me3_rep2"],
              row_cluster = False,
              col_cluster = False,
              row_colors=row_colors,
               robust=True,
               cmap = sns.cubehelix_palette(start=2.5, rot=0, dark=.25, light=.95, reverse=False, as_cmap=True),
              yticklabels  = False,
              figsize = (5,9),
               cbar_pos=(1, 0.65, 0.04, 0.18),
              cbar_kws={'ticks': [0.0, 3, 6]});

from matplotlib.patches import Patch
handles = [Patch(facecolor=lut[name]) for name in lut]

plt.legend(handles, lut, title='Clusters',
           bbox_to_anchor=(1.1, 0.62), bbox_transform=plt.gcf().transFigure, loc='upper right',frameon=False);

#Remove y label
ax = g.ax_heatmap
ax.set_ylabel("");

plt.savefig(f'{output_dir}ExtData_Fig1c_GBA_H3K4me3_differential_TSSs.pdf',
            bbox_inches="tight",
            transparent=True)
              


# #### Boxplots of H3K4me3 signals Across identified clusters

# In[ ]:


# Plot parameters
nrow = 2

ncol = 5

x_tickLabels = [string.split("_", maxsplit=1)[0] for string in exp_names ]

# Boxplot of H3K4me3 signal for each Replicate of ILC and for each cluster
fig, ax = plt.subplots(nrows=nrow,
                       ncols=ncol,
                      figsize = (16,6))

for df, cluster, ax in zip(cluster_df_long, cluster_names, ax.reshape(-1)):
    
    sns.boxplot(x="cellType", y="k4me3_signal", data=df, ax =ax)
    
    ax.set_xticklabels(x_tickLabels)
    
    ax.set_xlabel("")
    
    ax.set_title(f'Cluster {cluster}')
    

sns.despine()

fig.tight_layout()


# In[ ]:


#Isolate Cluster 3 and 4 that represent the ILC1 and ILC2 specific promoters
ilc1_ilc2_promoters = cluster_df_long[2:4]

#Remove ILC3 signal
ilc1_ilc2_promoters = [df[~df['cellType'].str.contains("ILC3")] for df in ilc1_ilc2_promoters]

#plot distribution of promoter-K4me3 signal 
nrow = 1

ncol = 2

x_tickLabels = [string.split("_", maxsplit=1)[0] for string in exp_names]

x_tickLabels = x_tickLabels[:4]

cluster_names = ["ILC1-specific \n H3K4me3 marked promoters",
              "ILC2-specific \n H3K4me3 marked promoters"]

fig, ax = plt.subplots(nrows=nrow,
                       ncols=ncol,
                      figsize = (8,4))

for df, cluster, ax in zip(ilc1_ilc2_promoters, cluster_names, ax.reshape(-1)):
    
      
    sns.boxplot(x="cellType",
                y="k4me3_signal",
                data=df,
                ax =ax,
                linewidth = 0.5)

    ax.set_xticklabels(x_tickLabels)

    ax.set_xlabel("")

    ax.set_title(f'Cluster {cluster}')
   
    sns.despine()

fig.tight_layout()

plt.savefig(f'{output_dir}ExtData_Fig1d_Fig1e_GBA_H3K4me3_differential_TSS_boxplots.pdf',
            bbox_inches="tight",
            transparent=True)


# #### Test if difference in H3K4me3 signal is statistically significant between ILC1 and ILC2

# In[ ]:


ilc1_cluster_data = pd.pivot(ilc1_ilc2_promoters[0], columns='cellType',values='k4me3_signal')


# In[ ]:


ilc2_cluster_data = pd.pivot(ilc1_ilc2_promoters[1], columns='cellType',values='k4me3_signal')


# In[ ]:


# ILC1_rep1 vs. ILC2_rep1 at promoters enriched with H3K4me3 in ILC1s
mannwhitneyu(x = ilc1_cluster_data['ILC1_GBA_H3K4me3_rep1'],
             y = ilc1_cluster_data['ILC2_GBA_H3K4me3_rep1'],
             alternative = 'greater',
            nan_policy = "omit")


# In[ ]:


# ILC1_rep2 vs. ILC2_rep2 at promoters enriched with H3K4me3 in ILC1s
mannwhitneyu(x = ilc1_cluster_data['ILC1_GBA_H3K4me3_rep2'],
             y = ilc1_cluster_data['ILC2_GBA_H3K4me3_rep2'],
             alternative = 'greater',
            nan_policy = "omit")


# In[ ]:


# ILC1_rep1 vs. ILC2_rep1 at promoters enriched with H3K4me3 in ILC2s
mannwhitneyu(x = ilc2_cluster_data['ILC1_GBA_H3K4me3_rep1'],
             y = ilc2_cluster_data['ILC2_GBA_H3K4me3_rep1'],
             alternative = 'less',
            nan_policy = "omit")


# In[ ]:


# ILC1_rep2 vs. ILC2_rep2 at promoters enriched with H3K4me3 in ILC2s
mannwhitneyu(x = ilc2_cluster_data['ILC1_GBA_H3K4me3_rep2'],
             y = ilc2_cluster_data['ILC2_GBA_H3K4me3_rep2'],
             alternative = 'less',
            nan_policy = "omit")


# #### Write BED files for ILC1 and ILC2 specific TSSs

# In[ ]:


#ILC1 signature BED file
cluster_df[2].reset_index().to_csv(f'{output_dir}ILC1_specific_k4me3-marked_TSSs.bed',
                                   sep = '\t',
                                   columns = ['level_0','level_1','level_2'],
                                  index = False,
                                  header = None)

#ILC2 signature BED file
cluster_df[3].reset_index().to_csv(f'{output_dir}ILC2_specific_k4me3-marked_TSSs.bed',
                                   sep = '\t',
                                   columns = ['level_0','level_1','level_2'],
                                  index = False,
                                  header = None)



