# Generates Extended Data figure 3i

library(ggrepel)
library(tidyverse)
library(Seurat)
library(Signac)


# Output Directory -------------------------------------------------------------
out.dir <- "./Figure_3/"


# Read Replicate 1 Data - Seurat Object ----------------------------------------
ilc.integrated <- readRDS(file = "./Figure_3/inVitro_ILC2_pcILC2_integrated_rep1.rds")

Idents(object = ilc.integrated) <- "seurat_clusters"

# Identify differential Motifs between cluster 9 (ILC1-like) and 
# ILC2s (clusters - 2,3,4,5,6,7)
clust.9.diff.motifs <- FindMarkers(ilc.integrated,
                                   ident.1 = 9,
                                   ident.2 = c(2,3,4,5,6,7),
                                   min.pct = 0.25,
                                   assay = "chromvar")

# Convert motif IDs to Motif names
clust.9.diff.motifs$motifName <- motif.name <- ConvertMotifID(ilc.integrated, id = row.names(clust.9.diff.motifs))


# Cluster 9 (ILC1-like) vs. ILC2 Motif Volcano Plot ----------------------------

clust9.diff.motif.plot <- clust.9.diff.motifs %>% 
  mutate(sig.color  = ifelse(-log10(p_val_adj) > -log10(0.05), '#1F9FE1', '#E0601E')) %>%
  arrange(sig.color)

extdata_3i <- ggplot(clust9.diff.motif.plot, aes(x=avg_log2FC, y=-log10(p_val_adj),  col = sig.color, label = motifName)) + 
  geom_point() +
  scale_color_identity() +
  geom_text_repel(box.padding = .2,
                  size=3,
                  segment.size = 0.25,
                  force_pull=1,
                  max.overlaps=10,
                  min.segment.length= 0.01,
                  data=subset(clust9.diff.motif.plot, (motifName == "TBX21" | motifName == "GATA3" | motifName == "RUNX3" |
                                                         motifName == "EOMES"| motifName == "STAT1" |
                                                         motifName == "SPIB"
                  )),
                  colour = 'black') +
  theme_classic() + 
  theme(legend.position="none") +
  ylab("-log10(p-value)") +
  xlab("Motif Enrichment \n log2(ILC1-like/ILC2)")

ggsave(extdata_3i,
       file = paste0(out.dir, "ExtData_Fig3i_ILC1-like_v_ILC2_diffMotif_volcanoPlot.pdf"),
       height = 4.5,
       width = 4.5,
       device = "pdf")

# Save Session Info ------------------------------------------------------------
writeLines(capture.output(sessionInfo()), "./sessionInfo/ExtData_Figure_3i_sessionInfo.txt")
