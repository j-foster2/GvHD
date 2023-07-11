# Generates Extended Data Figure 2k-l


library(Seurat)
library(Signac)
library(ggplot2)
library(tidyr)
library(dplyr)
library(stringr)


# Output Directory -------------------------------------------------------------
out.dir <- "./Figure_2/"


# Read Post-transplant Seurat Object -------------------------------------------
post_transplant <- readRDS(paste0(out.dir,"Post-transplant_SeuratObject_replicate1.rds"))


# Read Integrated Seurat Object ------------------------------------------------
ilc.combined <- readRDS(paste0(out.dir,"Post-transplant_replicate1_GBA_integration.rds"))


# Calculate Average ILC1/2/3 LDG signal on integration -------------------------
#ILC1 Lineage Defining genes
ilc1.genes <- c("sct_Tbx21","sct_Ifng","sct_Il21r","sct_Ccl5","sct_Ccl4",
                "sct_Ccl3","sct_Xcl1")

#ILC2 Lineage Defining genes
ilc2.genes <- c("sct_Gata3","sct_Hes1","sct_Lmo4","sct_Klf4","sct_Areg",
                "sct_Ccl1","sct_Csf2","sct_Il4", "sct_Il5", "sct_Il13")

#ILC3 Lineage Defining genes
ilc3.genes <- c("sct_Rorc","sct_Tcf7","sct_Foxs1","sct_Batf3","sct_Il22",
                "sct_Cx3cl1","sct_Il17f","sct_Gpx1")

innate.lymphoid.genes <- list(ilc1.genes, ilc2.genes, ilc3.genes)

ilc.type <- c("1", "2", "3")

for (j in 1:length(innate.lymphoid.genes)) {
  
  #Extract RNA-seq signal
  gene.signal = FetchData(object = ilc.combined, vars = innate.lymphoid.genes[[j]])
  #
  #Calculate per cell average of ILC1 genes
  gene.signal.avg = gene.signal %>%
    mutate(ILC.avg = rowMeans(.))
  
  #Add ILC average signal to meta.data data frame
  ilc.combined@meta.data[ , ncol(ilc.combined@meta.data) + 1] <- gene.signal.avg$ILC.avg
  
  #Update column name
  colnames(ilc.combined@meta.data)[ncol(ilc.combined@meta.data)] <- paste0("ILC", ilc.type[j], ".genes.avg")
  
}


# Annotate each cell as ILC1, 2, or, 3 -----------------------------------------

# Method #1 - Avg. Rank 

# Select normalized RNA signal for all lineage defining genes (across subtype) 
ilc.genes <- c(ilc1.genes, ilc2.genes, ilc3.genes)

ilc.genes.data <- FetchData(object = ilc.combined, vars = ilc.genes)

# Average Rank for each gene set
ilc.genes.rank <- as.data.frame(t(apply(ilc.genes.data, 1, rank)))

ilc.genes.rank$ILC1 <- rowMeans(ilc.genes.rank[,ilc1.genes])

ilc.genes.rank$ILC3 <- rowMeans(ilc.genes.rank[,ilc3.genes])

ilc.genes.rank$ILC2 <- rowMeans(ilc.genes.rank[,ilc2.genes])

ilc.genes.rank.avg <- ilc.genes.rank[,c("ILC1","ILC3","ILC2")]


# Annotate each cell with sub type label with highest rank
ilc.genes.rank.avg$cell.annno <-colnames(ilc.genes.rank.avg)[max.col(ilc.genes.rank.avg,ties.method="last")]

ilc.genes.rank.avg$ties <- (ilc.genes.rank.avg["ILC1"] == ilc.genes.rank.avg[["ILC2"]]) & (ilc.genes.rank.avg["ILC1"] == ilc.genes.rank.avg[["ILC3"]])

ilc.genes.rank.avg <- within(ilc.genes.rank.avg, cell.annno[ties == TRUE] <- 'Undetermined')

# Add data back to main SeuratObject
ilc.combined@meta.data$ILC.class <- ilc.genes.rank.avg$cell.annno


# Data Visualization -----------------------------------------------------------

#ILC Classification UMAP
ilc.classification.plot = DimPlot(ilc.combined,
                                  reduction = "umap",
                                  group.by="ILC.class",
                                  split.by = "orig.ident",
                                  order = c("ILC1","ILC2","ILC3","Undetermined"),
                                  combine = TRUE) + 
  scale_color_manual(values=c("#D3D3D3","#469990" ,"#892603", "#3A63AD")) +
  aes(stroke = 0.25) +
  ggtitle("")

ilc.classification.plot$data$orig.ident <- factor(ilc.classification.plot$data$orig.ident, levels = c("Gury_BenAri", "post_transplant"))

ggsave(ilc.classification.plot,
       file = paste0(out.dir, "ExtData_Fig2k_post-transpalnt_GBA_ILC_classification.pdf"),
       height = 4.6,
       width = 11,
       device = "pdf")


# Differential RNA and Motif Analysis ------------------------------------------

Idents(object = ilc.combined) <- "orig.ident"

# Subset Seurat Object for Post-transplant cells
ilc.combined.motif <- subset(x = ilc.combined, idents = "post_transplant")

# Transfer the integrated UMAP embedding to original post_transplant 
# Seurat Object
post_transplant@reductions <- ilc.combined.motif@reductions

# Transfer ILC classifications to original ex.ILC2 Seurat Object
post_transplant@meta.data[["ILC.class"]] <- ilc.combined.motif@meta.data[["ILC.class"]]

#Normalize RNA data
post_transplant <- SCTransform(post_transplant, assay = "RNA")

#Identify differential genes 
markers_rna <- presto:::wilcoxauc.Seurat(X = post_transplant,
                                         group_by = 'ILC.class',
                                         assay = 'data',
                                         seurat_assay = 'SCT',
                                         groups_use = c("ILC1","ILC2", "Undetrmined"))

#Identify differential regions of chromatin accessibility
markers_motifs <- presto:::wilcoxauc.Seurat(X = post_transplant,
                                            group_by = 'ILC.class',
                                            assay = 'data',
                                            seurat_assay = 'chromvar',
                                            groups_use = c("ILC1","ILC2", "Undetrmined"))
#Format differential data table
motif.names <- markers_motifs$feature

colnames(markers_rna) <- paste0("RNA.", colnames(markers_rna))

colnames(markers_motifs) <- paste0("motif.", colnames(markers_motifs))

markers_rna$gene <- markers_rna$RNA.feature

markers_motifs$gene <- ConvertMotifID(post_transplant, id = motif.names, assay = "ATAC")


celltype = "ILC1"

#Function to identify top expressed differential motif TF
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

#Table of candidate ILC1-like regulators
top.ILC1.motifs <- head(topTFs("ILC1"), 10)

top.ILC1.motifs.gene <- top.ILC1.motifs$gene

#Table of candidate ILC2 regulators
top.ILC2.motifs <- head(topTFs("ILC2"), 8)

top.ILC2.motifs.gene <- top.ILC2.motifs$gene


# Plot log2FC(motif) vs. log2FC(expression) ------------------------------------
topILC1.rna <- markers_rna %>% 
  dplyr::filter(RNA.group == "ILC1") %>%
  dplyr::filter(gene %in% c(top.ILC1.motifs.gene, top.ILC2.motifs.gene)) %>%
  arrange(gene)

topILC1.motif <- markers_motifs %>% 
  dplyr::filter(motif.group == "ILC1") %>%
  dplyr::filter(gene %in% c(toupper(top.ILC1.motifs.gene), toupper(top.ILC2.motifs.gene)) |
                  gene %in% c(top.ILC1.motifs.gene, top.ILC2.motifs.gene)) %>%
  arrange(gene)

topILC1.rna.motif <- cbind(topILC1.rna,topILC1.motif)

#Remove duplicate "gene" header 
names(topILC1.rna.motif)[length(names(topILC1.rna.motif))]<-"motif.gene" 

# Scatter Plot
Summary_ILC1_regulators <- topILC1.rna.motif %>% ggplot(aes(x=motif.logFC, y=RNA.logFC, size=RNA.avgExpr)) + 
  geom_point() +
  geom_text(label=topILC1.rna.motif$RNA.feature, position = position_nudge(y = -0.11), size = 5) +
  theme_classic() +
  geom_hline(yintercept=0, linetype="dashed", linewidth=0.5) +
  geom_vline(xintercept=0, linetype="dashed", linewidth=0.5) +
  xlab("Motif Enrichment \n log2(ILC1-like / ILC2)") + ylab("RNA Abundance \n log2(ILC1-like / ILC2)") +
  ggtitle("Candiate Regulators of ILC1-like Cells")

ggsave(Summary_ILC1_regulators,
       file = paste0(out.dir, "ExtData_Fig2l_candidate_regulator_scatter_summary.pdf"),
       height = 6,
       width = 8,
       device = "pdf")


# Save Session Info ------------------------------------------------------------
writeLines(capture.output(sessionInfo()), "./sessionInfo/ExtData_Figure_2k_l_sessionInfo.txt")

