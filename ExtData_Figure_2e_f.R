# Generate Extended Data Figure 2e-f


library(Seurat)
library(Signac)
library(ggplot2)
library(tidyr)
library(dplyr)
# library(stringr)


# Output Directory -------------------------------------------------------------
out.dir <- "./Figure_2/"


# Read ILC RNA Integrated Seurat Object ----------------------------------------
ilc.integrated <- readRDS(paste0(out.dir,"Pre_Post_transplant_multiome_integration.rds"))


# Visualize RNA integration Data -----------------------------------------------

# Integration UMAP
ilc.integrated$orig.ident <- factor(ilc.integrated$orig.ident,
                                    levels = c("pre_rep1","pre_rep2",'post_rep1','post_rep2'))

ilc.integrated.umap <- DimPlot(ilc.integrated,
                               reduction = "umap",
                               group.by = "orig.ident",
                               label = FALSE,
                               label.size = 2.5,
                               repel = TRUE,
                               cols = c("#005aa2",
                                        "#8772ce",
                                        "#88B570",
                                        "#238654")) + 
  aes(stroke = 0.25) +
  ggtitle("") +
  ylab("UMAP 2")+
  xlab("UMAP 1") +
  theme(legend.position = c(0.75, 0.20))

ggsave(ilc.integrated.umap,
       file = paste0(out.dir,"ExtData_Fig2e_pre_post_integration_allReps.pdf"),
       height = 4.5,
       width = 4.5,
       device = "pdf")



# Integration UMAP Split by sample
ilc.integrated$orig.ident <- factor(ilc.integrated$orig.ident,
                                    levels = c("pre_rep1","post_rep1",'pre_rep2','post_rep2'))

ilc.integrated.umap.split <- DimPlot(ilc.integrated,
                                     reduction = "umap",
                                     split.by = "orig.ident",
                                     group.by = "orig.ident",
                                     label = FALSE,
                                     label.size = 2.5,
                                     repel = TRUE,
                                     cols = c("#88B570",
                                              "#005aa2",
                                              "#61c5c4",
                                              "#8772ce",
                                              "#a77f48",
                                              "#C9575B",
                                              "#cf7231")) + aes(stroke = 0.25)

ggsave(ilc.integrated.umap.split,
       file = paste0(out.dir,"ExtData_Fig2e_pre_post_integration_allReps_split.pdf"),
       height = 5,
       width = 16,
       device = "pdf")


# Save Session Info ------------------------------------------------------------
writeLines(capture.output(sessionInfo()), "./sessionInfo/ExtData_Figure_2e_f_sessionInfo.txt")
