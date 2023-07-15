# Generates Figure 3i

library(tidyverse)
library(Seurat)
library(Signac)


# Output Directory -------------------------------------------------------------
out.dir <- "./Figure_3/"


# Read Replicate 1 Data - Seurat Object ----------------------------------------
ilc.integrated <- readRDS(file = "./Figure_3/inVitro_ILC2_pcILC2_integrated_rep1.rds")


# Figure 3i --------------------------------------------------------------------
Idents(object = ilc.integrated) <- "orig.ident"

cluster <- "pcILC2s"

# Regulator identification
markers_rna <- presto:::wilcoxauc.Seurat(X = ilc.integrated,
                                         group_by = 'orig.ident',
                                         assay = 'data',
                                         seurat_assay = 'SCT')

markers_motifs <- presto:::wilcoxauc.Seurat(X = ilc.integrated,
                                            group_by = 'orig.ident',
                                            assay = 'data',
                                            seurat_assay = 'chromvar')

motif.names <- markers_motifs$feature

colnames(markers_rna) <- paste0("RNA.", colnames(markers_rna))

colnames(markers_motifs) <- paste0("motif.", colnames(markers_motifs))

markers_rna$gene <- markers_rna$RNA.feature

markers_motifs$gene <- ConvertMotifID(ilc.integrated,
                                      id = motif.names, assay = "ATAC")

#a simple function to implement the procedure above
topTFs <- function(celltype, padj.cutoff = .05) {
  ctmarkers_rna <- dplyr::filter(
    markers_rna, RNA.group == celltype, RNA.padj < padj.cutoff, RNA.logFC > 0) %>% 
    arrange(-RNA.auc)
  ctmarkers_motif <- dplyr::filter(
    markers_motifs, motif.group == celltype, motif.padj < padj.cutoff, motif.logFC > 0) %>% 
    arrange(-motif.auc)
  
  ctmarkers_motif$gene <- tolower(ctmarkers_motif$gene)
  
  ctmarkers_motif$gene <- str_to_title(ctmarkers_motif$gene)
  
  top_tfs <- inner_join(
    x = ctmarkers_rna[, c(2, 11, 6, 7)], 
    y = ctmarkers_motif[, c(2, 1, 11, 6, 7)], by = "gene"
  )
  top_tfs$avg_auc <- (top_tfs$RNA.auc + top_tfs$motif.auc) / 2
  top_tfs <- arrange(top_tfs, -avg_auc)
  return(top_tfs)
}

clusterX.regulators <- head(topTFs(cluster),n = 10)

fig3i <- clusterX.regulators %>% 
  ggplot(aes(x = avg_auc,y = gene)) +
  geom_bar(stat = "identity")+
  scale_fill_gradient()+
  aes(y = fct_reorder(gene,avg_auc)) +
  theme_classic() +
  xlab("Average AUC") + 
  ylab("Candiate Regulators") +
  ggtitle("") +
  geom_text(aes(label= round(avg_auc, digits = 2)), vjust=1, hjust=1.25, color="white", size=3.5)


ggsave(fig3i,
       file =paste0(out.dir, "Fig3i_pcILC2_regulators.pdf"),
       height = 4,
       width = 4,
       device = "pdf")


# Save Session Info ------------------------------------------------------------
writeLines(capture.output(sessionInfo()), "./sessionInfo/Figure_3i_sessionInfo.txt")

