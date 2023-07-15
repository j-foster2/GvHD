# Generates Extended Data Figure 3h

library(ggrepel)
library(tidyverse)
library(Seurat)
library(Signac)


# Output Directory -------------------------------------------------------------
out.dir <- "./Figure_3/"


# Read Replicate 1 Data - Seurat Object ----------------------------------------
ilc.integrated <- readRDS(file = "./Figure_3/inVitro_ILC2_pcILC2_integrated_rep1.rds")


#Identify differential Genes between cluster 9 (ILC1-like) and ILC2s (clusters - 2,3,4,5,6,7)
clust.9.diff.genes <- FindMarkers(ilc.integrated,
                                  ident.1 = 9,
                                  ident.2 = c(2,3,4,5,6,7),
                                  min.pct = 0.25,
                                  assay = "SCT")

clust.9.diff.genes$geneName <- row.names(clust.9.diff.genes)


# Calculate BaseMean expression for differential genes -------------------------
clusters <- c(9,2,3,4,5,6,7)

cells.in.clusters <- ilc.integrated@meta.data[which(ilc.integrated@meta.data$seurat_clusters %in% clusters),] %>%
  row.names()

gene_test <- paste0("sct_",clust.9.diff.genes$geneName)

baseMean.data <-colMeans(FetchData(object = ilc.integrated,
                                   cells = cells.in.clusters,
                                   vars = gene_test))

clust.9.diff.genes$baseMean <- baseMean.data

clust.9.diff.genes$log2BaseMean <- log2(baseMean.data + 1)


# MA Plot ---------------------------------------------------------------------- 
clust.9.diff.genes <- clust.9.diff.genes %>% mutate(sig.color  = ifelse(-log10(p_val_adj) > -log10(0.05), '#1F9FE1', '#E0601E')) %>%
  arrange(sig.color)

extdata.fig3h <- ggplot(clust.9.diff.genes, aes(x=log2BaseMean, y=avg_log2FC,  col = sig.color, label = geneName)) +
  geom_point() +
  scale_color_identity() +
  geom_text_repel(box.padding = .25,
                  size=3,
                  segment.size = 0.25,
                  force_pull=1,
                  max.overlaps=30,
                  min.segment.length= 0.01,
                  data=subset(clust.9.diff.genes, (geneName == "Gata3" | geneName == "Il13" |
                                                     geneName == "Il5" | geneName == "Ifng" |
                                                     geneName == "Gzmb" | geneName == "Inpp4b"
                  )),
                  colour = 'black') +
  theme_classic() +
  theme(legend.position="none") +
  ylab("-log10(p-value)") +
  xlab("RNA Abundance \n log2(ILC1-like / ILC2)")

ggsave(extdata.fig3h,
       file = paste0(out.dir, "ExtData_Fig3h_cluster9_v_ILC2s_MAPlot.pdf"),
       height = 4.5,
       width = 4.5,
       device = "pdf")

# Save Session Info ------------------------------------------------------------
writeLines(capture.output(sessionInfo()), "./sessionInfo/ExtData_Figure_3h_sessionInfo.txt")

