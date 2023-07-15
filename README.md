# Multiomic Single Cell Evaluation Reveals Inflammatory Cytokines Affect Innate Lymphoid Cell Fate After Allogeneic Stem Cell Transplantation

## Directories 

(Code needs to be exected as described below for all directories to appear)

- `processed_data`: data run through data-appropriate pipeline (i.e. CellRanger)
- `CellRanger`: example code of how cellranger was run on single nucleus multiome and scATAC datasets  
- `sessionInfo`: sessionInfo for each of the R scripts
- `Bruce_Processed_Data`: processed data from Bruce et al. (PMID: 28375154)
- `GBA_Processed_Data`: MARS-seq processed data from Gury-BenAri et al. (PMID: 27545347)
- `referenceData`: mm10 gene coordinates
- `enrichrData`: output from Enrichr analyses 
- `jupyterNotebooks`: Store notebook versions of scripts used to processed data and generate figures
- `Figure_1`: output from scripts for figure 1
- `Figure_2`: output from scripts for figure 2
- `Figure_3`: output from scripts for figure 3
- `Figure_4`: output from scripts for figure 4

## Data

Create a folder called `processed_data`, and download the processed files from GEO SuperSeries `GSE232003`.

## Figure 1 workflow

- [generate_geneLevel_RNA_expression_Fig1b.ipynb](generate_geneLevel_RNA_expression_Fig1b.ipynb) - generates input file to convert Bruce et al transcript-level RNA estimates  to gene-level estimates
- [transcritToGeneConversion_mm10.R](transcritToGeneConversion_mm10.R) - generates gene-level RNA Abunance estimates for Bruce et al RNA abundance data
- [bruce_RNA_laurie_foster_K4me3_data.ipynb](bruce_RNA_laurie_foster_K4me3_data.ipynb) - generates dataframe rows all genes columns include average  gene-level RNA abundance (Bruce et al.) and average H3K4me3 signal proximal to TSSs
- [Figure_1b.R](Figure_1b.R) - generates figure 1b plot
- [generate_GBA_H3K4me3_count_matrix.R](generate_GBA_H3K4me3_count_matrix.R) - generates H3K4me3 count matrix at TSSs (including 300 bp upstream and 500 downstream)  
- [GBA_H3K4me3_DESeq2.R](GBA_H3K4me3_DESeq2.R) - identifes TSSs with differential H3K4me3 signal across three subtypes of ILCs
- [ExtData_Figure_1.ipynb](ExtData_Figure_1.ipynb) - hierarchical clustering on differential H3K4me3 singal and generates Extended Data Figure 1c-e 
- [Figure_1c.ipynb](Figure_1c.ipynb) - generates figure 1c plot

## Figure 2 workflow

- [integrate_multiome_replicate1.R](integrate_multiome_replicate1.R) - integrates first replicate of mouse single-nucleus multiome data
- [integrate_multiome_replicate2.R](integrate_multiome_replicate2.R) - integrates second replicate of mouse single-nucleus multiome data
- [Figure_2b.R](Figure_2b.R) - generates figure 2b plot 
- [Figure_2c.R](Figure_2c.R) - generates figure 2c plot
- [Figure_2d.R](Figure_2d.R) - generates figure 2d plot
- [ExtData_Figure_2d.R](ExtData_Figure_2d.R) - generates extended data figure 2d plot
- [post_transplant_geneSet_Enrichr.ipynb](post_transplant_geneSet_Enrichr.ipynb) - description of post-transplant gene set enrichr analysis and data formatting
- [Figure_2e.R](Figure_2e.R) - generates figure 2e plot
- [Figure_2f.R](Figure_2f.R) - generates figure 2f plot
- [Figure_2g.R](Figure_2g.R) - generates figure 2g plot
- [Figure_2h.R](Figure_2h.R) - generates figure 2h plot
- [Figure_2i_m.R](Figure_2i_m.R) - generates figure 2i-m plots
- [integrate_post_transplant_replicate2_GBA.R](integrate_post_transplant_replicate2_GBA.R) - integrates second replicate of post-transplant multiome data with GBA MARS-seq data
- [Figure_2n_o.R](Figure_2n_o.R) - generates figure 2n-o plots
- [integrate_all_multiome_data.R](integrate_all_multiome_data.R) - integrates all single-nucleus multiome data using RNA signal
- [ExtData_Figure_2e_f.R](ExtData_Figure_2e_f.R) - generates extended data figure 2e-f plots
- [ExtData_Figure_2g_j.R](ExtData_Figure_2g_j.R) - generates extended data figure 2g-j plots
- [integrate_post_transplant_replicate1_GBA.R](integrate_post_transplant_replicate1_GBA.R) - integrates first replicate of post-transplant multiome data with GBA MARS-seq data
- [ExtData_Figure_2k_l.R](ExtData_Figure_2k_l.R) - generates extended data figure 2k_l plots

## Figure 3 workflow

- [Figure_3d.R](Figure_3d.R) - generates figure 3d plot
- [mouse_ATAC_peak_union_set.py](mouse_ATAC_peak_union_set.py) - identifies the union set of ATAC peaks (Figure 3e)
- [generate_pcILC2_ILC2_ATAC_count_matrix.R](generate_pcILC2_ILC2_ATAC_count_matrix.R) - counts of all peaks in union set (Figure 3e)
- [ILC2_pcILC2_ATAC_DESeq2.R](ILC2_pcILC2_ATAC_DESeq2.R) - differential analysis of peaks in union set (Figure 3e)
- [Figure_3e.py](Figure_3e.py) - generates figure 3e plot
- [integrate_inVitro_ILC2_pcILC2_replicate1.R](integrate_inVitro_ILC2_pcILC2_replicate1.R) - integrates first replicate of mouse ILC2 and pcILC2 data
- [Figure_3g.R](Figure_3g.R) - generates figure 3g plot
- [Figure_3h.R](Figure_3h.R) - generates figure 3h plot
- [Figure_3i.R](Figure_3i.R) - generates figure 3i plot
- [Figure_3j.R](Figure_3j.R) - generates figure 3j plot
- [integrate_inVitro_ILC2_pcILC2_replicate2.R](integrate_inVitro_ILC2_pcILC2_replicate2.R) - integrates second replicate of mouse ILC2 and pcILC2 data

## Figure 4 workflow
