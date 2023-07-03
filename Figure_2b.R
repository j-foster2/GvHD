# Creates the plot shown in figure 2b. The UMAP reflects the dimension
# reduction performed on the integrated data.This script also calculates the
# percentage of pre-transplant cells that cluster with the post-transplant cells
# and vice versa. 

library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)

# Output Directory -------------------------------------------------------------

out.dir <- "./Figure_2/"


# Read Integrated Seurat Object ------------------------------------------------

ilc.integrated <- readRDS(paste0(out.dir,"pre-transplant_post-transplant_integration_replicate2.rds"))


# FIG. 2b - UMAP annotated with sample ID -------------------------------------

# Set order of samples
ilc.integrated$orig.ident <- factor(ilc.integrated$orig.ident, levels = c("pre_transplant","post_transplant"))

ilc.integrated.umap <- DimPlot(ilc.integrated,
                               reduction = "umap",
                               group.by = "orig.ident",
                               label = FALSE,
                               label.size = 2.5,
                               repel = TRUE,
                               cols = c("#2C5D73", "#7cd37b")) + 
  aes(stroke = 0.25) +
  ggtitle("") + 
  theme(legend.position = c(0.01, 0.15))

ggsave(ilc.integrated.umap,
       file = paste0(out.dir,"Fig2b_pre-post_transplant_integration_umap.pdf"),
       height = 4.5,
       width = 4.5 ,
       device = "pdf")

# Counting Numbers of cells in FIG. 2b -----------------------------------------

#Number of post-transplant cells that are in pre-transplant clusters
num_post_in_pre <- ilc.integrated@meta.data %>% 
  dplyr::filter(orig.ident == "post_transplant" &
                  (integrated_snn_res.0.4 == 0 | integrated_snn_res.0.4 == 1 |integrated_snn_res.0.4 == 2)) %>% 
  nrow()

num_post <- ilc.integrated@meta.data %>% 
  dplyr::filter(orig.ident == "post_transplant" ) %>% nrow()

perct_post_in_pre <- num_post_in_pre / num_post

writeLines(capture.output(perct_post_in_pre), "./Figure_2/Fig2b_num_postCells_cluster_with_preCells.txt")

#Number of pre-transplant cells that are in post-transplant clusters
num_pre_in_post <- ilc.integrated@meta.data %>% 
  dplyr::filter(orig.ident == "pre_transplant" &
                  (integrated_snn_res.0.4 == 3 | integrated_snn_res.0.4 == 4 |integrated_snn_res.0.4 == 5)) %>% 
  nrow()

num_pre <- ilc.integrated@meta.data %>% 
  dplyr::filter(orig.ident == "pre_transplant" ) %>% nrow()

perct_pre_in_post <- num_pre_in_post / num_pre

writeLines(capture.output(perct_pre_in_post), "./Figure_2/Fig2b_num_preCells_cluster_with_postCells.txt")
