# Generates Figure 2n-o


library(Seurat)
library(Signac)
library(ggplot2)
library(tidyverse)
library(stringr)


# Output Directory -------------------------------------------------------------
out.dir <- "./Figure_2/"


# Read Post-transplant Seurat Object -------------------------------------------
ilc.integrated <- readRDS(paste0(out.dir,"pre-transplant_post-transplant_integration_replicate1.rds"))


# Assign cluster names ---------------------------------------------------------

# Use UAMP and cluster overlay to assign labels to clusters
ilc.integrated.umap <- DimPlot(ilc.integrated,
                               reduction = "umap",
                               group.by = "integrated_snn_res.0.3",
                               split.by = "orig.ident",
                               label = TRUE,
                               label.size = 2.5,
                               repel = TRUE) + 
  aes(stroke = 0.25) +
  ggtitle("") + 
  theme(legend.position = c(0.01, 0.15))


# Label data
Idents(object = ilc.integrated) <- "integrated_snn_res.0.3"

ilc.integrated <- RenameIdents(object = ilc.integrated,
                               `0` = "Pre-1",
                               `1` = "Pre-3",
                               `2` = "Post-1",
                               `3` = "Pre-2",
                               `4` = "Pre-4",
                               `5` = "Pre-5")


# Extended data Fig2g ----------------------------------------------------------

ilc.integrated.umap <- DimPlot(ilc.integrated,
                               reduction = "umap",
                               group.by = "orig.ident",
                               label = FALSE,
                               label.size = 2.5,
                               repel = TRUE,
                               cols = c("#2C5D73", "#7cd37b")) + 
  aes(stroke = 0.25) +
  ggtitle("") + 
  theme(legend.position = c(0.65, 0.90))

ggsave(ilc.integrated.umap,
       file = paste0(out.dir,"ExtData_Fig2g_pre_post_transplant_integration.pdf"),
       height = 4.5,
       width = 4.5,
       device = "pdf")


# Extended data Fig2h ----------------------------------------------------------

ilc.integrated$orig.ident <- factor(ilc.integrated$orig.ident,
                                    levels = c("pre_transplant","post_transplant"))

Idents(ilc.integrated) <- factor(x = Idents(ilc.integrated), levels = c("Pre-1", "Pre-2","Pre-3",
                                                                        "Pre-4", "Pre-5", "Post-1"))

ilc.integrated.umap.split <- DimPlot(ilc.integrated,
                                     reduction = "umap",
                                     split.by = "orig.ident",
                                     label = F,
                                     label.size = 2.5,
                                     repel = TRUE,
                                     cols = c("#005aa2",
                                              "#c7c657",
                                              "#62c6c5",
                                              "#8372b3",
                                              "#ce722f",
                                              "#88B570")) + aes(stroke = 0.25) +
  ggtitle("") 

ggsave(ilc.integrated.umap.split,
       file = paste0(out.dir,"ExtData_Fig2h_pre_post_transplant_integration.pdf"),
       height = 4,
       width = 8,
       device = "pdf")



# Extended Data Figi-j ---------------------------------------------------------

#Notes: 
# TF candidate analysis where comparison is only on the pre- or post-transplant
# groups.

# Select clustering resolution
resolution <- "integrated_snn_res.0.3"

Idents(object = ilc.integrated) <- resolution

# Store cluster names to iterate through
clusters.names <- levels(ilc.integrated)

# Store number of motifs at gained regions of chromatin accessibility
num.motifs.at.gained <- list()

# Store number of putative regulators identified for each cluster 
num.putative.regulators <- list()

ilc.integrated <- RenameIdents(object = ilc.integrated,
                               `0` = "Pre-1",
                               `1` = "Pre-3",
                               `2` = "Post-1",
                               `3` = "Pre-2",
                               `4` = "Pre-4",
                               `5` = "Pre-5")

cluster.labels <- c("Pre-1", "Pre-3", "Post-1",
                    "Pre-2", "Pre-4", "Pre-5")

colors <- c("#005aa2", "#62c6c5", "#88B570",
            "#c7c657", "#8372b3", "#ce722f")

for (i in 1:length(clusters.names)) {
  cluster = clusters.names[i]
  
  ##### Putative regulators of pre-transplant clusters #####
  if (cluster == "0" | cluster == "1" | cluster == "3" | cluster == "4" | cluster == "5") {
    
    # Differential RNA analysis 
    markers_rna <- presto:::wilcoxauc.Seurat(X = ilc.integrated,
                                             group_by = resolution,
                                             assay = 'data',
                                             seurat_assay = 'SCT',
                                             groups_use = c(cluster, "2"))
    
    # Differential motif analysis 
    markers_motifs <- presto:::wilcoxauc.Seurat(X = ilc.integrated,
                                                group_by = resolution,
                                                assay = 'data',
                                                seurat_assay = 'chromvar',
                                                groups_use = c(cluster, "2"))
    
    # Store motif names
    motif.names <- markers_motifs$feature
    
    # Annotate differential RNA data with "RNA. "prefix - track RNA information 
    colnames(markers_rna) <- paste0("RNA.", colnames(markers_rna))
    
    # Annotate differential motif data with "motif. "prefix - track motif information
    colnames(markers_motifs) <- paste0("motif.", colnames(markers_motifs))
    
    # Add column with gene names
    markers_rna$gene <- markers_rna$RNA.feature
    
    # Add column with understandable motif name 
    markers_motifs$gene <- ConvertMotifID(ilc.integrated, id = motif.names, assay = "ATAC")
    
    # Function to identify candidate regulators 
    # TF with increased RNA and whose motif is present
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
    
    # Top 10 putative regulator for a given cluster
    clusterX.regulators <- head(topTFs(cluster), n = 10)
    
    # Store number of putative regulator for a given cluster
    num.putative.regulators <- append(num.putative.regulators, list(nrow(topTFs(cluster))))
    
    # Bar plot of top 10 putative regulator for a given cluster
    clusterX.regulators %>% 
      ggplot(aes(x = avg_auc,y = gene)) +
      geom_bar(stat = "identity", fill = colors[i])+
      scale_fill_gradient()+
      aes(y = fct_reorder(gene,avg_auc)) +
      theme_classic() +
      xlab("Average AUC") + 
      ylab("Candiate Regulators") +
      ggtitle(paste0("Cluster ",cluster)) + 
      scale_x_continuous(limits = c(0,1), breaks= c(0, 0.5, 1), labels=c(0, 0.5, 1)) +
      geom_text(aes(label= round(avg_auc, digits = 2)), vjust=1.25, hjust=1.25, color="white", size=3.5)
    
    output.filename = paste0(out.dir,"ExtData_Fig2i_", cluster.labels[i], "_cluster", cluster, "_v_postClusters_candidate_TF.pdf")
    
    ggsave(file = output.filename,
           height = 4.6,
           width = 3,
           device = "pdf")
  }
  
  ##### Putative regulators of post-transplant clusters #####
  else if (cluster == "2") {
    
    # Differential RNA analysis 
    markers_rna <- presto:::wilcoxauc.Seurat(X = ilc.integrated,
                                             group_by = resolution,
                                             assay = 'data',
                                             seurat_assay = 'SCT',
                                             groups_use = c(cluster, "0","1","3","4","5"))
    
    # Differential motif analysis 
    markers_motifs <- presto:::wilcoxauc.Seurat(X = ilc.integrated,
                                                group_by = resolution,
                                                assay = 'data',
                                                seurat_assay = 'chromvar',
                                                groups_use = c(cluster, "0","1","3","4","5"))
    
    # Store motif names
    motif.names <- markers_motifs$feature
    
    # Annotate differential RNA data with "RNA. "prefix - track RNA information 
    colnames(markers_rna) <- paste0("RNA.", colnames(markers_rna))
    
    # Annotate differential motif data with "motif. "prefix - track motif information
    colnames(markers_motifs) <- paste0("motif.", colnames(markers_motifs))
    
    # Add column with gene names
    markers_rna$gene <- markers_rna$RNA.feature
    
    # Add column with understandable motif name 
    markers_motifs$gene <- ConvertMotifID(ilc.integrated, id = motif.names, assay = "ATAC")
    
    # Function to identify candidate regulators 
    # TF with increased RNA and whose motif is present
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
    
    # Top 10 putative regulator for a given cluster
    clusterX.regulators <- head(topTFs(cluster),n = 10)
    
    # Store number of putative regulator for a given cluster
    num.putative.regulators <- append(num.putative.regulators, list(nrow(topTFs(cluster))))
    
    # Bar plot of top 10 putative regulator for a given cluster
    clusterX.regulators %>% 
      ggplot(aes(x = avg_auc,y = gene)) +
      geom_bar(stat = "identity", fill = colors[i])+
      scale_fill_gradient()+
      aes(y = fct_reorder(gene,avg_auc)) +
      theme_classic() +
      xlab("Average AUC") + 
      ylab("Candiate Regulators") +
      ggtitle(paste0("Cluster ",cluster)) + 
      scale_x_continuous(limits = c(0,1), breaks= c(0, 0.5, 1), labels=c(0, 0.5, 1)) +
      geom_text(aes(label= round(avg_auc, digits = 2)), vjust=1.25, hjust=1.25, color="white", size=3.5)
    
    
    output.filename = paste0(out.dir,"ExtData_Fig2j_", cluster.labels[i], "_cluster", cluster, "_v_preClusters_candidate_TF.pdf")
    
    ggsave(file = output.filename,
           height = 4.6,
           width = 3,
           device = "pdf")  
  }
}

# Save Session Info ------------------------------------------------------------

writeLines(capture.output(sessionInfo()), "./sessionInfo/ExtData_Figure_2g_j_sessionInfo.txt")

