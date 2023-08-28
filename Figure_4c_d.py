#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# Generates Figures 4c-d


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


# In[ ]:


out_dir = "./Figure_4/"


# In[ ]:


hILC2_peakFile = f'{out_dir}hILC2-specific_ATAC_peaks.bed'

pc_hILC2_peakFile = f'{out_dir}pc-hILC2-specific_ATAC_peaks.bed'

background_peakFile = f'{out_dir}hILC2_pc-hILC2_unchanged_ATAC_peaks.bed'


# #### ATAC signal at GATA3 

# In[ ]:


# GATA3 motif file 
gata3_motif_file = "./Figure_4/hILC2-specific_peak_motifs/knownResults/known23.motif"


# In[ ]:


get_ipython().run_cell_magic('bash', '-s "$out_dir" "$hILC2_peakFile" "$gata3_motif_file"', '\nannotatePeaks.pl \\\n$2 \\\nhg38 \\\n-m $3 > $1hILC2-specific_peak_motifs/knownResults/ATAC_peaks_GATA3_motif.txt\n')


# In[ ]:


# Read hILC2 annotated ATAC peak bed file 
hILC2_motifs = pd.read_csv("./Figure_4/hILC2-specific_peak_motifs/knownResults/ATAC_peaks_GATA3_motif.txt",
                         sep = '\t')


# In[ ]:


# Filter for ATAC peaks with GATA3 motif
gata3_coord = hILC2_motifs[hILC2_motifs['GATA3(Zf)/iTreg-Gata3-ChIP-Seq(GSE20898)/Homer Distance From Peak(sequence,strand,conservation)'].notnull()]

# Filter for genomic coordinates of ATAC peaks with GATA3 motif
gata3_coord = gata3_coord[['Chr', 'Start', 'End']]

# Write coordinates to bed file

gata3_coord.to_csv(f"{out_dir}hILC2_ATAC_Peaks_with_GATA3_motif.bed",
                  sep = '\t',
                  index = False,
                  header = None)



# In[ ]:


# Read hILC2 and pc-hILC2 bigWig files
hILC2_bw = sorted(glob.glob('./processed_data/hILC2_ATAC_rep*openChromatin.bw'))

pc_hILC2_bw = sorted(glob.glob('./processed_data/pc-hILC2_ATAC_rep*openChromatin.bw'))

fn_atac_bw = hILC2_bw + pc_hILC2_bw

#Convert bigWig file name list to space delimited string
str_fn_atac_bw = " ".join(fn_atac_bw)

#Gather experiment names 
exp_names = ['hILC2_rep1','hILC2_rep2',
            'pc-hILC2_rep1','pc-hILC2_rep2']

#Convert experiment names list to space delimited string
str_exp_names = " ".join(exp_names)


# In[ ]:


# Store BED file of ATAC peaks with GATA3 motif

bedWithGata3 = f"{out_dir}hILC2_ATAC_Peaks_with_GATA3_motif.bed"


# In[ ]:


get_ipython().run_cell_magic('bash', '-s  "$str_fn_atac_bw" "$str_exp_names" "$bedWithGata3" "$out_dir"', '\n#Define Variables\nbw_files=$1\nexp_names=$2\nbed_file=$3\noutDir=$4\n\ncomputeMatrix reference-point \\\n--scoreFileName $bw_files \\\n--regionsFileName $bed_file \\\n--beforeRegionStartLength 1000 \\\n--afterRegionStartLength 1000 \\\n--referencePoint center \\\n--samplesLabel $exp_names \\\n--outFileName $outDir\\hILC2_peaks_with_GATA3Peaks.npz \\\n--outFileNameMatrix $outDir\\hILC2_peaks_with_GATA3Peaks.txt \\\n-p 8\n')


# In[ ]:


# Parameters
window_size = 100

data_division = len(exp_names)


# Read normalized ATAC signal at peaks with GATA3 motif
atacSignalWithMotif = np.loadtxt(f'{out_dir}hILC2_peaks_with_GATA3Peaks.txt',
                                      dtype='float64',
                                      skiprows=3)

# Split matrix based on number of samples
atacSignalWithMotifSplit = np.hsplit(atacSignalWithMotif, data_division)

# Columnwise Mean to create average line plot
atacSignalMeanGATA3 = [] 

#used to center the data around appropriate window
bp_range = list(range(-window_size, window_size)) 

for i in range(len(atacSignalWithMotifSplit)):
    atacSignalMeanGATA3.append(atacSignalWithMotifSplit[i].mean(axis=0))
    
    
fix, ax = plt.subplots(figsize=(5,5))

ax.plot(bp_range,atacSignalMeanGATA3[0], label = 'hILC2_rep1', color = "#395e73")
ax.plot(bp_range,atacSignalMeanGATA3[1], label = 'hILC2_rep1', color = "#395e73", linestyle = "--")
ax.plot(bp_range,atacSignalMeanGATA3[2], label = 'pc-hILC2_rep1', color = "#e39a44")
ax.plot(bp_range,atacSignalMeanGATA3[3], label = 'pc-hILC2_rep2', color = "#e39a44", linestyle = "--")

ax.set_ylabel("Normalized ATAC Signal")

ax.legend(loc=1, ncol=1, borderaxespad=0)

sns.despine()

plt.setp((ax), xticks=[-window_size,0,window_size], 
         xticklabels=['-{}kb'.format(int(window_size/100)),'center','{}kb'.format(int(window_size/100))])


# #### ATAC signal at Tbet (pc-hILC2) 

# In[ ]:


# Tbet motif file 
Tbet_motif_file = "./Figure_4/pc_hILC2-specific_peak_motifs/homerResults/motif5.motif"



# In[ ]:


get_ipython().run_cell_magic('bash', '-s "$out_dir" "$pc_hILC2_peakFile" "$Tbet_motif_file"', '\nannotatePeaks.pl \\\n$2 \\\nhg38 \\\n-m $3 > $1pc_hILC2-specific_peak_motifs/homerResults/ATAC_peaks_Tbet_motif.txt\\\n')


# In[ ]:


# Read exILC2 annotated ATAC peak bed file 
Tbet_motifs = pd.read_csv("./Figure_4/pc_hILC2-specific_peak_motifs/homerResults/ATAC_peaks_Tbet_motif.txt",
                         sep = '\t')


# In[ ]:


# Filter for ATAC peaks with GATA3 motif
Tbet_coord = Tbet_motifs[Tbet_motifs['7-YNTMACACCT,BestGuess:Tbet(T-box)/CD8-Tbet-ChIP-Seq(GSE33802)/Homer(0.952) Distance From Peak(sequence,strand,conservation)'].notnull()]

# Filter for genomic coordinates of ATAC peaks with GATA3 motif
Tbet_coord = Tbet_coord[['Chr', 'Start', 'End']]

# Write coordinates to bed file

Tbet_coord.to_csv(f"{out_dir}pc-hILC2_ATAC_Peaks_with_Tbet_motif.bed",
                  sep = '\t',
                  index = False,
                  header = None)



# In[ ]:


# Read hILC2 and pc-hILC2 bigWig files
hILC2_bw = sorted(glob.glob('./processed_data/hILC2_ATAC_rep*openChromatin.bw'))

pc_hILC2_bw = sorted(glob.glob('./processed_data/pc-hILC2_ATAC_rep*openChromatin.bw'))

fn_atac_bw = hILC2_bw + pc_hILC2_bw

#Convert bigWig file name list to space delimited string
str_fn_atac_bw = " ".join(fn_atac_bw)

#Gather experiment names 
exp_names = ['hILC2_rep1','hILC2_rep2',
            'pc-hILC2_rep1','pc-hILC2_rep2']

#Convert experiment names list to space delimited string
str_exp_names = " ".join(exp_names)


# In[ ]:


# Store BED file of ATAC peaks with Tbet motif

bedWithTbet = f"{out_dir}pc-hILC2_ATAC_Peaks_with_Tbet_motif.bed"


# In[ ]:


bedWithTbet


# In[ ]:


get_ipython().run_cell_magic('bash', '-s  "$str_fn_atac_bw" "$str_exp_names" "$bedWithTbet" "$out_dir"', '\n#Define Variables\nbw_files=$1\nexp_names=$2\nbed_file=$3\noutDir=$4\n\ncomputeMatrix reference-point \\\n--scoreFileName $bw_files \\\n--regionsFileName $bed_file \\\n--beforeRegionStartLength 1000 \\\n--afterRegionStartLength 1000 \\\n--referencePoint center \\\n--samplesLabel $exp_names \\\n--outFileName $outDir\\pc_hILC2_peaks_with_TbetPeaks.npz \\\n--outFileNameMatrix $outDir\\pc_hILC2_peaks_with_TbetPeaks.txt \\\n-p 8\n')


# In[ ]:


# Parameters
window_size = 100

data_division = len(exp_names)

# Read normalized ATAC signal at peaks with GATA3 motif
atacSignalWithMotif = np.loadtxt(f'{out_dir}pc_hILC2_peaks_with_TbetPeaks.txt',
                                      dtype='float64',
                                      skiprows=3)

# Split matrix based on number of samples
atacSignalWithMotifSplit = np.hsplit(atacSignalWithMotif, data_division)

# Columnwise Mean to create average line plot
atacSignalMean = [] 

#used to center the data around appropriate window
bp_range = list(range(-window_size, window_size)) 

for i in range(len(atacSignalWithMotifSplit)):
    atacSignalMean.append(atacSignalWithMotifSplit[i].mean(axis=0))
    
    
fix, ax = plt.subplots(figsize=(5,5))

ax.plot(bp_range,atacSignalMean[0], label = 'ILC2_rep1', color = "#395e73")
ax.plot(bp_range,atacSignalMean[1], label = 'ILC2_rep1', color = "#395e73", linestyle = "dotted")
ax.plot(bp_range,atacSignalMean[2], label = 'exILC2_rep1', color = "#e39a44")
ax.plot(bp_range,atacSignalMean[3], label = 'exILC2_rep2', color = "#e39a44", linestyle = "dotted")

ax.set_ylabel("Normalized ATAC Signal")

ax.legend(loc=1, ncol=1, borderaxespad=0)

sns.despine()

plt.setp((ax), xticks=[-window_size,0,window_size], 
         xticklabels=['-{}kb'.format(int(window_size/100)),'center','{}kb'.format(int(window_size/100))])


# In[ ]:


#
fix, ax = plt.subplots(1,2,figsize=(10,5), sharey=True)

ax[0].plot(bp_range,atacSignalMeanGATA3[0], label = 'ILC2_rep1', color = "#395e73")
ax[0].plot(bp_range,atacSignalMeanGATA3[1], label = 'ILC2_rep1', color = "#395e73", linestyle = "dotted")
ax[0].plot(bp_range,atacSignalMeanGATA3[2], label = 'exILC2_rep1', color = "#e39a44")
ax[0].plot(bp_range,atacSignalMeanGATA3[3], label = 'exILC2_rep2', color = "#e39a44", linestyle = "dotted")


ax[0].set_ylabel("Normalized ATAC Signal")

plt.setp((ax[0]), xticks=[-window_size,0,window_size], 
         xticklabels=['-{}kb'.format(int(window_size/100)),'center','{}kb'.format(int(window_size/100))])

ax[1].plot(bp_range,atacSignalMean[0], label = 'ILC2_rep1', color = "#395e73")
ax[1].plot(bp_range,atacSignalMean[1], label = 'ILC2_rep1', color = "#395e73", linestyle = "dotted")
ax[1].plot(bp_range,atacSignalMean[2], label = 'exILC2_rep1', color = "#e39a44")
ax[1].plot(bp_range,atacSignalMean[3], label = 'exILC2_rep2', color = "#e39a44", linestyle = "dotted")


ax[1].legend(loc=1, ncol=1, borderaxespad=0)

sns.despine()

plt.setp((ax[1]), xticks=[-window_size,0,window_size], 
         xticklabels=['-{}kb'.format(int(window_size/100)),'center','{}kb'.format(int(window_size/100))])


plt.savefig(f'{out_dir}Fig4d_hILC2_pc-hILC2_ATAC_Peaks_with_GATA3_TBET_motif_avgLinePlot.pdf'
            ,bbox_inches="tight",
            transparent=True)


# #### Visualization of Homer de novo results

# In[ ]:


#Paths to deNovo Homer Output 

hILC2_path = f'{out_dir}hILC2-specific_peak_motifs/homerResults/'

pc_hILC2_path = f'{out_dir}pc_hILC2-specific_peak_motifs/homerResults/'


# In[ ]:


def homerBarPlot(path, numMotifs):
      
    tf_names = []
    pvalues = []
    
    for num in range(1, numMotifs + 1):
        
        filename = f'motif{num}.motif'
        
        with open (path + filename) as file:
            first_line = file.readline().split('\t')
            
        #Isolate transcriptoin factor
        tf_names.append(first_line[1].split('/')[0].split(":")[-1])
    
        #Isolate p-value 
        pval = first_line[-1].split(',')[-1].rstrip('\n').lstrip('P:').split("-")[1]
    
        pvalues.append(int(pval))
    
    #Format data for Seaborn barplot
    dict_motif_pval = {'Motif': tf_names,'pval':pvalues}

    df_motif_pval = pd.DataFrame(data = dict_motif_pval)
    
    #Create Barplot of Figure
    fig, ax = plt.subplots(figsize=(2.5, 4))

    sns.barplot(x = 'pval',
                y = 'Motif',
                data = dict_motif_pval,
                palette = "Blues_r",
                ax = ax)
    
    ax.set_xlabel('Motifs Enrichment \n -log$_{10}$ p-value', 
                 fontweight = "bold")
    
    sns.despine()
    
    return(fig)


# In[ ]:


#ILC2 de novo barplot 
ILC2_barplot = homerBarPlot(hILC2_path, 5)

ILC2_barplot.savefig(f'{out_dir}Fig4c_hILC2_homer_deNovo_barPlot.pdf',            
            bbox_inches="tight",
            transparent=True
           )


# In[ ]:


#exILC2 de novo barplot 
pc_hILC2_barplot = homerBarPlot(pc_hILC2_path, 5)

pc_hILC2_barplot.savefig(f'{out_dir}Fig4c_pc-hILC2_homer_deNovo_barPlot.pdf',            
            bbox_inches="tight",
            transparent=True
           )

