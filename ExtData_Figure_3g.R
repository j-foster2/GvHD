library(tidyverse)
library(Seurat)
library(Signac)
library(stringr)
library(ggplot2)


# Output Directory -------------------------------------------------------------
out.dir <- "./Figure_3/"


# Read Replicate 1 Data - Seurat Object ----------------------------------------
ilc.integrated <- readRDS(file = "./Figure_3/inVitro_ILC2_pcILC2_integrated_rep1.rds")


# Extended Data Figure 3g ------------------------------------------------------

#RNA Integration UMAP
ilc.integrated.umap <- DimPlot(ilc.integrated,
                               reduction = "umap",
                               group.by = "seurat_clusters",
                               label = TRUE,
                               label.size = 2.5,
                               repel = TRUE,
                               cols = c("#88B570",
                                        "#005aa2",
                                        "#61c5c4",
                                        "#8772ce",
                                        "#a77f48",
                                        "#C9575B",
                                        "#cf7231",
                                        "#43bf7b",
                                        "#bc541b",
                                        "#cebf32")) + 
  aes(stroke = 0.25) +
  ggtitle("")

ggsave(ilc.integrated.umap, file = paste0(out.dir,"ExtData_Fig3g_ILC2_pcILC2_UMAP_seuratClusters.pdf"),
       height = 4.5, width = 5 , device = "pdf")

# Save Session Info ------------------------------------------------------------
writeLines(capture.output(sessionInfo()), "./sessionInfo/ExtData_Figure_3g_sessionInfo.txt")

