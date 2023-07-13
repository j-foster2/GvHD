# Generates Figure 3g.

library(Seurat)
library(Signac)
library(tidyverse)


# Output Directory -------------------------------------------------------------
out.dir <- "./Figure_3/"


# Read Replicate 1 Data - Seurat Object ----------------------------------------
ilc.integrated <- readRDS(file = "./Figure_3/inVitro_ILC2_pcILC2_integrated_rep1.rds")


# Figure 3g --------------------------------------------------------------------
ilc.integrated.umap <- DimPlot(ilc.integrated,
                               reduction = "umap",
                               group.by = "orig.ident",
                               label = FALSE,
                               label.size = 2.5,
                               repel = TRUE,
                               cols = c("#2C5D73", "#F0972B")) + 
  aes(stroke = 0) +
  ggtitle("") +
  theme(legend.position = c(0.01, 0.15))

ggsave(ilc.integrated.umap,
       file = paste0(out.dir,"Fig3g_ILC2_pcILC2_UMAP.pdf"),
       height = 4.5,
       width = 4.5,
       device = "pdf")


# Save Session Info ------------------------------------------------------------
writeLines(capture.output(sessionInfo()), "./sessionInfo/Figure_3g_sessionInfo.txt")

