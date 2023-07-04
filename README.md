# Multiomic Single Cell Evaluation Reveals Inflammatory Cytokines Affect Innate Lymphoid Cell Fate After Allogeneic Stem Cell Transplantation

## Directories 

(Code needs to be exected as described below for all directories to appear)

- `processed_data`: run through data-appropriate pipeline (i.e. CellRanger)
- `CellRanger`: example code of hour cellranger was run on single nucleus multiome and scATAC datasets  
- `sessionInfo`: sessionInfo for each of the R scripts
- `Figure_1`: output from scripts for figure 1
- `Figure_2`: output from scripts for figure 2
- `Figure_3`: output from scripts for figure 3
- `Figure_4`: output from scripts for figure 4


## Data

Create a folder called `processed_data`, and download the processed files from GEO SuperSeries `GSE232003`.

## Figure 1 workflow

- [generate_geneLevel_RNA_expression_Fig1b.ipynb](generate_geneLevel_RNA_expression_Fig1b.ipynb) - generate input file to convert Bruce et al transcript-level RNA estimates  to gene-level estimates
- [transcritToGeneConversion_mm10.R](transcritToGeneConversion_mm10.R) - generate gene-level RNA Abunance estimates for Bruce et al RNA abundance data
- [bruce_RNA_laurie_foster_K4me3_data.ipynb](bruce_RNA_laurie_foster_K4me3_data.ipynb) - generate dataframe rows all genes columns include average  gene-level RNA abundance (Bruce et al.) and average H3K4me3 signal proximal to TSSs
- [Figure_1b.R](Figure_1b.R) - generate figure 1b plot

## Figure 2 workflow

- [integrate_multiome_replicate1.R](integrate_multiome_replicate1.R) - integrate first replicate of mouse single-nucleus multiome data
- [integrate_multiome_replicate2.R](integrate_multiome_replicate2.R) - integrate second replicate of mouse single-nucleus multiome data
- [Figure_2b.R](Figure_2b.R) - generate figure 2b plot 
- [Figure_2c.R](Figure_2c.R) - generate figure 2c plot

## Figure 3 workflow

## Figure 4 workflow
