#!/usr/bin/env python
# coding: utf-8

# ### Notes
# #### Generation of Enrichr Data
# 1. Post-transplant Data Key
#     1. ./Figure_2/ILC2_transplant_differential_genes_cluster_3.txt --> Post1
#     2. ./Figure_2/ILC2_transplant_differential_genes_cluster_4.txt --> Post2
#     3. ./Figure_2/ILC2_transplant_differential_genes_cluster_5.txt --> Post3
#     
# 2. Gene lists were passed to Enrichr via online portal (https://maayanlab.cloud/Enrichr/)
# 
# 3. Terms called using Mouse Gene Atlas Database were input for this script
#  
#  
# #### Description of Script  
# 1. Pulled the top5 terms for each of the post-transplant gene sets
# 
# 2. Calculated overlap, if the term was not in the top five of a given gene set then overlap was set to 0

# In[ ]:


import pandas as pd
import numpy as np
import glob
import os


# In[ ]:


# Output Directory
output_dir = "./Figure_2/"


# In[ ]:


# read Enrichr results (Mouse Gene Atlas)
enrichr_files = glob.glob("./enrichrData/Mouse_Gene_Atlas_table_Clust_Post*txt")


# In[ ]:


# Isolate cluster names from file name 
cluster_Names = [os.path.basename(file).split("_")[-1].split(".")[0] for file in enrichr_files]


# In[ ]:


enrichr_data = [pd.read_csv(file, sep ='\t') for file in enrichr_files]

# Annotate each dataset with cluster
enrichr_data = [df.assign(cluster = cluster) for df,cluster in zip(enrichr_data,cluster_Names)]

# Isolate the top 5 terms
enrichr_data = [df.head(5) for df in enrichr_data]

# Select Term, adj P-value, cluster
enrichr_data = [df[['Term','Overlap','Adjusted P-value','cluster']] for df in enrichr_data]

# Calculate overlap percentage
overlap_pct = [df['Overlap'].apply(lambda x: int(x.split("/")[0]) / int(x.split("/")[1])) for df in enrichr_data]

for df, overlap in zip(enrichr_data, overlap_pct):
    df['Overlap'] = overlap 

enrichr_data  = pd.concat(enrichr_data, ignore_index=True)

# Copied for appending overlap to final matrix
overlap_data = enrichr_data


# In[ ]:


# Long-to-Wide
enrichr_data = overlap_data.pivot(index='cluster', columns='Term', values='Adjusted P-value')

# nan to 1
enrichr_data = enrichr_data.fillna(1)

# -log10 transform data
enrichr_data = -np.log10(enrichr_data)

# replace -0 with 0 values
enrichr_data = enrichr_data.replace(-0.0, 0) 


# In[ ]:


# wide to long
enrichr_data = pd.melt(enrichr_data.reset_index(), id_vars='cluster',value_vars=enrichr_data.columns.to_list())


# In[ ]:


# Add overlap pct to the matrix
enrichr_data = enrichr_data.assign(overlap = float(0))

for idx, row in overlap_data.iterrows():
    for idx2, row2 in enrichr_data.iterrows():
        if row['Term'] == row2['Term'] and row['cluster'] == row2['cluster']:
            enrichr_data.at[idx2,'overlap'] = row['Overlap']


# In[ ]:


# Save Matrix
enrichr_data.to_csv(f"{output_dir}post_transplant_enrichr.txt",
                   sep = '\t',
                   index = False)

