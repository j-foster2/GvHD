#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# Calculates the H3K4me3 signal at every TSS across genome and associates with  
# RNA abunance of the corresponding genes. 


# In[ ]:


import glob
import os
import pandas as pd


# #### Output Directory

# In[ ]:


# Store path to output directory
output_dir = "./Figure_1/"


# #### Create BED file of Refseq promotor coordinates

# In[ ]:


#Read mm10 RefSeq file 
mm10_refseq = pd.read_csv('./referenceData/mm10_refseq_long_transcript_nohaplo_2020-01-25.bed',
                          delimiter='\t',
                          header=None)
#Format RefSeq file
mm10_refseq.columns = ['chr','start','end','strand','geneName','transcriptID']

mm10_refseq.drop(columns=['transcriptID'], inplace = True)

#Create TSS BED file (-300,+500)
upstream_dist = 300
downstream_dist = 500

with open(f'{output_dir}/refseq_promotor.bed', 'w') as output:
    for idx, row in mm10_refseq.iterrows():
        if row[3] == "+" and row[3] != 'strand':
            output.write('%s\t%s\t%s\t%s\n' % \
                         (row[0], row[1] - upstream_dist, row[1] + downstream_dist,row[4]))
        elif row[3] == "-" and row[3] != 'strand':
            output.write('%s\t%s\t%s\t%s\n' % \
                         (row[0], row[2] - downstream_dist,row[2] + upstream_dist,row[4]))


# #### Read gene-level RNA abundance estimates from Bruce et al. ()

# In[ ]:


# Read in gene-level RNA abundance data for Bruce ILC2s
rna_data = pd.read_csv('./Figure_1/ILC2_Bruce_geneLevel_RNA_abundance.txt',
                   sep = '\t')

rna_data = rna_data.rename(columns={"Unnamed: 0": "geneName"})

rna_data = rna_data.set_index("geneName")

rna_data = rna_data.loc[:,"abundance.ILC2_RNA_Bruce_rep1":"abundance.ILC2_RNA_Bruce_rep3"]


# #### Calculate K4me3 read-depth normalized signal at gene promoters

# In[ ]:


str_k4me3_bw


# In[ ]:


# Read the bigWig file names for all acute ChIP-seq data
path_To_Laurie_Foster_K4me3_Data = './processed_data/'

k4me3_bw = sorted(glob.glob(f'{path_To_Laurie_Foster_K4me3_Data}blfH3K4ME3_rep*'))

# Convert bigWig file name list to space delimited string
str_k4me3_bw = " ".join(k4me3_bw)

# Experiment Names
exp_names = ["H3K4me3_rep1", "H3K4me3_rep2"]

# Convert experiment names list to space delimited string
str_exp_names = " ".join(exp_names)


# In[ ]:


# Save promoter bed file 
mm10_promoter_bed = f'{output_dir}refseq_promotor.bed'


# In[ ]:


get_ipython().run_cell_magic('bash', '', '\n# Check version of deeptools\ndeeptools --version\n')


# In[ ]:


get_ipython().run_cell_magic('bash', '-s  "$str_k4me3_bw" "$str_exp_names" "$mm10_promoter_bed" "$output_dir"', '\nmultiBigwigSummary BED-file \\\n--bwfiles $1 \\\n--BED $3 \\\n--labels $2 \\\n--outFileName $4ILC2_H3K4me3_promoter.npz \\\n--outRawCounts $4ILC2_H3K4me3_promoter.tab \\\n-p 24\n')


# #### Join RNA and H3K4me3 DataFrames

# In[ ]:


# Annotate deepTools matrix with geneNames
k4me3_signal = pd.read_csv(f'{output_dir}ILC2_H3K4me3_promoter.tab',
                           sep = '\t',
                           header = 0,
                           names = ['chr','start','end','H3K4me3_rep1','H3K4me3_rep2'],
                           index_col = ['chr','start','end'])


# In[ ]:


# Join RNA and H3K4me3 data into one dataframe
mm10_promoter_bed_df = pd.read_csv(mm10_promoter_bed,
                                   sep = '\t',
                                   header = None,
                                   names = ['chr','start','end', 'geneName'],
                                   index_col = ['chr','start','end'])

mm10_promoter_bed_df = k4me3_signal.join(mm10_promoter_bed_df)

mm10_promoter_bed_df = mm10_promoter_bed_df.set_index('geneName')

rna_data_k4me3_data = rna_data.join(mm10_promoter_bed_df)

rna_data_k4me3_data = rna_data_k4me3_data.assign(ILC2_RNA_avg = rna_data_k4me3_data.loc[:,"abundance.ILC2_RNA_Bruce_rep1":"abundance.ILC2_RNA_Bruce_rep3"].mean(axis = 1))

rna_data_k4me3_data = rna_data_k4me3_data.assign(ILC2_K4me3_avg = rna_data_k4me3_data.loc[:,"H3K4me3_rep1":"H3K4me3_rep2"].mean(axis = 1))


# #### Save Data for Visualization in R

# In[ ]:


#Save dataframe for visualization in R
rna_data_k4me3_data.to_csv(f'{output_dir}bruce_rna_laurie_k4me3_data.txt',
                           sep = '\t',
                           index = False)


# In[ ]:




