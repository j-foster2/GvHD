#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# hILC2 homer analysis


# In[ ]:


import pandas as pd
import glob
import numpy as np
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
import seaborn as sns


# #### Generate Bed files for HOMER motif anlaysis

# In[ ]:


out_dir = "./Figure_4/"


# In[ ]:


ILC_clust_data = pd.read_csv(f"{out_dir}Fig4b_ILC_clust_data.txt",
                            sep = '\t',
                            index_col= "peakID")


# In[ ]:


ILC2_peaks = ILC_clust_data[ILC_clust_data['log2FoldChange'] > 0]

ILC2_peaks.to_csv(f'{out_dir}hILC2-specific_ATAC_peaks.bed',
                 sep = '\t',
                 columns = ['chr','start','end'],
                 header = None,
                 index = False)

pc_hILC2_peaks = ILC_clust_data[ILC_clust_data['log2FoldChange'] < 0]

pc_hILC2_peaks.to_csv(f'{out_dir}pc-hILC2-specific_ATAC_peaks.bed',
                 sep = '\t',
                 columns = ['chr','start','end'],
                 header = None,
                 index = False)


# In[ ]:


ILC_deseq_unchanged = pd.read_csv(f'{out_dir}hILC2_pc-hILC2_deseq2_results.txt',
                        sep = '\t',
                       index_col = "Unnamed: 0")

ILC_deseq_unchanged = ILC_deseq_unchanged[ILC_deseq_unchanged['padj'] >= 0.05]

ILC_union_bed = pd.read_csv(f'{out_dir}pc-hILC2_hILC2_union_200bp.bed',
                            sep = '\t',
                           index_col = 'peakName')

ILC_deseq_unchanged.join(ILC_union_bed).to_csv(f'{out_dir}hILC2_pc-hILC2_unchanged_ATAC_peaks.bed',
                                     sep = '\t',
                                     header = None,
                                     index = False,
                                     columns = ['Chr', 'Start','End'])


# In[ ]:


hILC2_peakFile = f'{out_dir}hILC2-specific_ATAC_peaks.bed'

pc_hILC2_peakFile = f'{out_dir}pc-hILC2-specific_ATAC_peaks.bed'

background_peakFile = f'{out_dir}hILC2_pc-hILC2_unchanged_ATAC_peaks.bed'


# #### Motifs at open ATAC peaks in hILC2s

# In[ ]:


get_ipython().run_cell_magic('bash', '-s "$out_dir" "$hILC2_peakFile" "$background_peakFile"', '\nfindMotifsGenome.pl \\\n$2 \\\nhg38 \\\n$1hILC2-specific_peak_motifs \\\n-size 200 \\\n-bg $3 \\\n-preparse \\\n-preparsedDir $1preparse_hg38\n')

