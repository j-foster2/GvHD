#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# Create a text file that associates experiment name and path to salmon file


# In[ ]:


import pandas as pd
import glob
import os
import datetime


# ### Output Directory 

# In[ ]:


output_dir = './Figure_1/'


# ### Collect  sf files

# In[ ]:


# Read Bruce et al. transcript level RNA abundance data 
pathToBruceData = './Bruce_Processed_Data/'

bruce_rna_sf = sorted(glob.glob(f'{pathToBruceData}*sf'))

# Create text file for conversation of transcript to gene-level RNA abundance values
with open(f'{output_dir}bruce_sf_files.txt', 'w') as output:
    for file_path in bruce_rna_sf:
        output.write("%s\n" % file_path)

# Experiment names 
exp_names = [os.path.basename(name) for name in bruce_rna_sf]

exp_names = [x.replace('.sf', '') for x in exp_names]


# ### Format input file for transcript to gene level count data

# In[ ]:


ilc_sf_files = pd.read_csv(f'{output_dir}bruce_sf_files.txt', sep = '\t',header = None)

ilc_sf_files['sample_names'] = exp_names

ilc_sf_files=ilc_sf_files.rename(columns = {0:'file_paths'})

ilc_sf_files.to_csv(f'{output_dir}bruce_sf_files_formatted.txt',sep = '\t',header = None, index = False, columns = ['sample_names','file_paths'])

