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
out_dir = "./Figure_3/"


# In[ ]:


# Create a BED file for GREAT analysis

ILC_deseq = pd.read_csv(f'{out_dir}pcILC2_ILC2_deseq2_results_diff_0.05.txt',
                        sep = '\t',
                        index_col = "Unnamed: 0")

ILC_union_bed = pd.read_csv(f'{out_dir}pcILC2_ILC2_union_200bp.bed',
                            sep = '\t',
                           index_col = 'peakName')

ILC_deseq.join(ILC_union_bed).to_csv(f'{out_dir}pcILC2_ILC2_differential_ATAC_peaks.bed',
                                     sep = '\t',
                                     header = None,
                                     index = False,
                                     columns = ['Chr', 'Start','End'])

ILC_deseq_GREAT = ILC_deseq.assign(peakID = [f'peak{i}' for i in range(ILC_deseq.shape[0])])

ILC_deseq_GREAT.join(ILC_union_bed).to_csv(f'{out_dir}pcILC2_ILC2_differential_ATAC_peaks_GREAT.bed',
                                     sep = '\t',
                                     header = None,
                                     index = False,
                                     columns = ['Chr', 'Start','End','peakID'])


# In[ ]:


pcILC2_bw = sorted(glob.glob('./processed_data/mILC1_ATAC_rep*bw'))

ILC2_bw = sorted(glob.glob('./processed_data/mILC2_ATAC_rep*bw'))

fn_chip_bw = pcILC2_bw + ILC2_bw

#Convert bigWig file name list to space delimited string
str_fn_chip_bw = " ".join(fn_chip_bw)

#Gather experiment names 
exp_names = ["ILC1_rep1", "ILC1_rep2",
            "ILC2_rep1", "ILC2_rep2", "ILC2_rep3"]

#Convert experiment names list to space delimited string
str_exp_names = " ".join(exp_names)


# In[ ]:


#ILC1 specific promoters 
ILC_diff_peak_bed = f'{out_dir}pcILC2_ILC2_differential_ATAC_peaks.bed'


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
# 1. pcILC2_ILC2_differential_ATAC_peaks.bed was run through GREAT online portal 
#     1. Parameters:
#         1. Species Assembly: Mouse: GRCm38 (UCSC mm10, Dec. 2011)
#         2. Associating genomic regions with genes: Single nearest Gene within 1000 kb
#         3. Include curated regulatory domains = True

# In[ ]:


# Annotate the deepTools matrix with proximal genes 

file = './greatData/pcILC2_ILC2_diff_ATAC_GREAT_geneAnno.txt'

ilc_diff_peaks = pd.read_csv(file, sep = '\t')

ilc_diff_peaks.columns = ['peakNAME', 'gene_dist', 'association_rule']

ilc_diff_peaks.drop(columns=['association_rule'], inplace=True)

# Drop peaks that have no associated gene
ilc_diff_peaks = ilc_diff_peaks[~ilc_diff_peaks['gene_dist'].str.contains('NONE')]

ilc_diff_peaks[['geneName', 'dist']] = ilc_diff_peaks['gene_dist'].str.split(" ", n = 1,expand = True)

ilc_diff_peaks.drop(columns=['gene_dist'], inplace=True)

# Handles situation where there are multiple genes for a given ROE
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

# # Save BED of differential regions annotated with proximal genes
# ILC_clust_data.to_csv(f"{out_dir}ILC2_pcILC2_differential_ATAC_regions_geneNames.bed",
#                      columns = ['chr','start','end','geneName'],
#                      index = False,
#                      sep = '\t')


# In[ ]:


#Create Matrix with proximal gene set as index
ILC_clust_matrix = ILC_clust_data.reset_index().set_index("geneName").loc[:,"ILC1_rep1":"ILC2_rep3"]


# In[ ]:


#Reorder the Matrix (for plotting)
ILC_clust_matrix = ILC_clust_matrix[['ILC2_rep1','ILC2_rep2','ILC2_rep3','ILC1_rep1','ILC1_rep2']]

# Update Column Names
# Note ILC1 was renamed to pcILC2 to more accurately reflect the nature 
# of cells. These cells are ILC2s treated with a cytokine cocktail. The
# pc prefix = proinflammatory cytokine 
ILC_clust_matrix = ILC_clust_matrix.rename(columns={'ILC1_rep1': 'pcILC2_rep1', 'ILC1_rep2': 'pcILC2_rep2'})


# #### Generate Plot

# In[ ]:


ilc_clustmap = sns.clustermap(ILC_clust_matrix,
              z_score=0, method="average",
              cmap = sns.diverging_palette(230, 50, as_cmap=True),
              col_cluster = False,
              yticklabels=False,
              figsize = (5,9),
              cbar_pos=(1., 0.1, 0.04, 0.18),
              cbar_kws={'ticks': [-1.5, 0.0, 1.5]})

plt.savefig(f'{out_dir}Fig3e_ILC2_pcILC2_differential_ATAC_clusterMap.pdf',
            bbox_inches="tight",
            transparent=True)


# In[ ]:




