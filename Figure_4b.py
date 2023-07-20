#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pandas as pd
import glob
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
import seaborn as sns


# #### Output Directory

# In[ ]:


#Create Output Dir variable
out_dir = "./Figure_4/"


# In[ ]:


# Create a BED file for GREAT analysis

ILC_deseq = pd.read_csv(f'{out_dir}hILC2_pc-hILC2_deseq2_results_diff_0.05.txt',
                        sep = '\t',
                        index_col = "Unnamed: 0")

ILC_union_bed = pd.read_csv(f'{out_dir}pc-hILC2_hILC2_union_200bp.bed',
                            sep = '\t',
                           index_col = 'peakName')

ILC_deseq.join(ILC_union_bed).to_csv(f'{out_dir}pc-hILC2_hILC2_differential_ATAC_peaks.bed',
                                     sep = '\t',
                                     header = None,
                                     index = False,
                                     columns = ['Chr', 'Start','End'])

ILC_deseq_GREAT = ILC_deseq.assign(peakID = [f'peak{i}' for i in range(ILC_deseq.shape[0])])

ILC_deseq_GREAT.join(ILC_union_bed).to_csv(f'{out_dir}pc-hILC2_hILC2_differential_ATAC_peaks_GREAT.bed',
                                     sep = '\t',
                                     header = None,
                                     index = False,
                                     columns = ['Chr', 'Start','End','peakID'])


# In[ ]:


diff_peaks = ILC_deseq.join(ILC_union_bed)
#Create ILC2 signature 
ilc2_sig = diff_peaks[diff_peaks['log2FoldChange'] > 0]

ilc2_sig = ilc2_sig.assign(strand = "*")

ilc2_sig.to_csv(f"{out_dir}huILC2_ATAC_sig.bed",
               sep = '\t',
               header = None,
               index = False,
               columns = ['Chr','Start','End','strand'])

#Create exILC2 signature
exilc2_sig = diff_peaks[diff_peaks['log2FoldChange'] < 0]

exilc2_sig = exilc2_sig.assign(strand = "*")

exilc2_sig.to_csv(f"{out_dir}huexILC2_ATAC_sig.bed",
               sep = '\t',
               header = None,
               index = False,
               columns = ['Chr','Start','End','strand'])



# In[ ]:


ILC2_bw = sorted(glob.glob('./processed_data/hILC2_ATAC_rep*bw'))

pcILC2_bw = sorted(glob.glob('./processed_data/pc-hILC2_ATAC_rep*bw'))

fn_chip_bw =   ILC2_bw + pcILC2_bw

#Convert bigWig file name list to space delimited string
str_fn_chip_bw = " ".join(fn_chip_bw)

#Gather experiment names 
exp_names = ["hILC2_rep1", "hILC2_rep2",
            "pc-hILC2_rep1", "pc-hILC2_rep2"]

#Convert experiment names list to space delimited string
str_exp_names = " ".join(exp_names)


# In[ ]:


# hILC2/pc-hILC2 union set
ILC_diff_peak_bed = f'{out_dir}pc-hILC2_hILC2_differential_ATAC_peaks.bed'


# In[ ]:


get_ipython().run_cell_magic('bash', '', '\ndeeptools --version\n')


# In[ ]:


get_ipython().run_cell_magic('bash', '-s  "$str_fn_chip_bw" "$str_exp_names" "$ILC_diff_peak_bed" "$out_dir"', '\n\nmultiBigwigSummary BED-file \\\n--bwfiles $1 \\\n--BED $3 \\\n--labels $2 \\\n--outFileName $4pcILC2_ILC2_diff_peaks_openChromatin_Signal.npz \\\n--outRawCounts $4pcILC2_ILC2_diff_peaks_openChromatin_Signal.tab \\\n-p 4\n')


# ### Hierarchical Clustering of mILC2 mExILC2 Signature ATAC Peaks

# In[ ]:


ILC_peak_signal = pd.read_csv(f'{out_dir}pcILC2_ILC2_diff_peaks_openChromatin_Signal.tab',
                                  sep = '\t',
                                  names = ['chr','start','end'] + exp_names,
                                 header = 0)


# #### Annotate ATAC peaks with proximal genes (GREAT)
# 1. hILC2_pc-hILC2_GREAT_basalPlusExtension.txt was run through GREAT online portal 
#     1. Parameters:
#         1. Species Assembly: Human: GRCh38 (UCSC hg38, Dec. 2013)
#         2. Associating genomic regions with genes: Basal plus extension 
#             1. Proximal 5.0 kb
#             2. Upstream 1.0 kb
#             3. Distal 1000 kb
#         3. Include curated regulatory domains = True
#         

# In[ ]:


# Annotate the deepTools matrix with proximal genes 

file = './greatData/hILC2_pc-hILC2_GREAT_basalPlusExtension.txt'

ilc_diff_peaks = pd.read_csv(file, sep = '\t')

ilc_diff_peaks.columns = ['peakNAME', 'gene_dist', 'association_rule']

ilc_diff_peaks.drop(columns=['association_rule'], inplace=True)

#Book Keeping
print('Number ATAC Peaks not associated with gene: {}'.
      format(ilc_diff_peaks[ilc_diff_peaks['gene_dist'].str.contains('NONE')].shape[0]))

#Drop ROE that have no associated gene
ilc_diff_peaks = ilc_diff_peaks[~ilc_diff_peaks['gene_dist'].str.contains('NONE')]

ilc_diff_peaks[['geneName', 'dist']] = ilc_diff_peaks['gene_dist'].str.split(" ", n = 1,expand = True)

ilc_diff_peaks.drop(columns=['gene_dist'], inplace=True)

#Handles situation where there are multiple genes for a given ROE
if ilc_diff_peaks.loc[ilc_diff_peaks['dist'].str.contains(","),'dist'].empty != True: 
    ilc_diff_peaks.loc[ilc_diff_peaks['dist'].str.contains(","),'dist'] =  ilc_diff_peaks.loc[ilc_diff_peaks['dist'].str.contains(","),'dist'].tolist()[0].split(",")[0]
    
    #Isolate the distance value to nearest gene
    ilc_diff_peaks['dist'] = ilc_diff_peaks['dist'].str.lstrip('(+(-').str.rstrip(' )')

    #Type cast Distance to nearest gene as int64
    ilc_diff_peaks['dist'] = ilc_diff_peaks['dist'].astype('int64')
else:   
   #Isolate the distance value to nearest gene
    ilc_diff_peaks['dist'] = ilc_diff_peaks['dist'].str.lstrip('(+(-').str.rstrip(' )')

    #Type cast Distance to nearest gene as int64
    ilc_diff_peaks['dist'] = ilc_diff_peaks['dist'].astype('int64')



# In[ ]:


# Annotate ATAC signal matrix with proximal gene Names
ILC_peak_singal = ILC_peak_signal.set_index(["chr","start","end"])

# Add coordinates to the DESeq2 results (ILC2_)
ILC_deseq_GREAT = ILC_deseq_GREAT.join(ILC_union_bed).reset_index().rename(columns = {"Chr":"chr", "Start":"start", "End":"end"}).set_index(["chr","start","end"])

#Add unique peak IDs to dataFrame
ILC_clust_data = ILC_peak_singal.join(ILC_deseq_GREAT).reset_index().set_index("peakID")

#Add gene names to df
ILC_clust_data = ILC_clust_data.join(ilc_diff_peaks.set_index("peakNAME"))



# In[ ]:


#Create Matrix with proximal gene set as index
ILC_clust_matrix = ILC_clust_data.reset_index().set_index("geneName").loc[:,"hILC2_rep1":"pc-hILC2_rep2"]


# In[ ]:


#Reorder the Matrix (for plotting)
ILC_clust_matrix = ILC_clust_matrix[["hILC2_rep1", "hILC2_rep2","pc-hILC2_rep1", "pc-hILC2_rep2"]]


# #### Generate Plot

# In[ ]:


ilc_clustmap = sns.clustermap(ILC_clust_matrix,
              z_score=0, method="average",
              cmap = sns.diverging_palette(230, 50, as_cmap=True),
              col_cluster = False,
              yticklabels=False,
              figsize = (5,9),
              cbar_pos=(1., 0.1, 0.04, 0.18))

plt.savefig(f'{out_dir}Fig4b_ILC2_pcILC2_differential_ATAC_clusterMap.pdf',
            bbox_inches="tight",
            transparent=True)


# In[ ]:


# Save data frame used for heatmap
ILC_clust_data.to_csv(f"{out_dir}Fig4b_ILC_clust_data.txt",
                     sep = '\t')

