# Generates plot for Figure 2c


library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)


# Output Directory -------------------------------------------------------------

out.dir <- "./Figure_2/"


# Read Integrated Seurat Object ------------------------------------------------

ilc.integrated <- readRDS(paste0(out.dir,"pre-transplant_post-transplant_integration_replicate2.rds"))


# FIG. 2c - UMAP -clust. res.0.4 - pre v. post clustering - Split by sample-----

# Cluster data overlayed on UMAP was derived from an analysis where the 
#resolution parameter was set to 0.4.

# Set order of samples
ilc.integrated$orig.ident <- factor(ilc.integrated$orig.ident, levels = c("pre_transplant","post_transplant"))

# Create UMAP
ilc.integrated.umap.split <- DimPlot(ilc.integrated,
                                     reduction = "umap",
                                     group.by = 'integrated_snn_res.0.4',
                                     split.by = "orig.ident",
                                     label = FALSE,
                                     label.size = 2.5,
                                     repel = TRUE,
                                     cols = c("#a77f48",
                                              "#005aa2",
                                              "#62c6c5",
                                              "#238654",
                                              "#88B570",
                                              "#56813e")) + aes(stroke = 0.25) +
  ggtitle("") 

# Save plot
ggsave(ilc.integrated.umap.split, file = paste0(out.dir,"Fig2c_pre-post_transplant_integration_umap_split_clustered.pdf"),
       height = 4, width = 8, device = "pdf")

# Save Session Info ------------------------------------------------------------

writeLines(capture.output(sessionInfo()), "./sessionInfo/Figure_2c_sessionInfo.txt")
