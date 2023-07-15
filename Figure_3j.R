# Generates Figure 3j


library(tidyverse)
library(Seurat)
library(Signac)


# Output Directory -------------------------------------------------------------
out.dir <- "./Figure_3/"


# Read Replicate 1 Data - Seurat Object ----------------------------------------
ilc.integrated <- readRDS(file = "./Figure_3/inVitro_ILC2_pcILC2_integrated_rep1.rds")


# Figure 3f --------------------------------------------------------------------

# Expression overlay for a couple of ICL1 associated genes
genesOfInterest <- c('Tbx21', 'Ifng')

for (gene in genesOfInterest) {
  
  
  gene_plot <- FeaturePlot(ilc.integrated,
                           features = paste0("sct_",gene),
                           reduction = 'umap',
                           order = T)+ 
    theme(legend.position = "bottom")
  
  motif.name <- ConvertMotifID(ilc.integrated, name = toupper(gene))
  
  ggsave(gene_plot,
         file = paste0(out.dir,"Fig3j_ILC2_pcILC2_",gene,"_RNA_expression_UMAP.pdf"),
         height = 4,
         width = 4,
         device = "pdf")
}


# Save Session Info ------------------------------------------------------------
writeLines(capture.output(sessionInfo()), "./sessionInfo/Figure_3j_sessionInfo.txt")
