# Generates plot for Figure 2I-M


library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)
library(chromVAR)   
library(JASPAR2020) #Works
library(TFBSTools)
library(motifmatchr)
library(BSgenome.Mmusculus.UCSC.mm10)
library(stringr)
library(forcats)
library(ggrepel)



# Output Directory -------------------------------------------------------------

out.dir <- "./Figure_2/"


# Read Integrated Seurat Object ------------------------------------------------

ilc.integrated <- readRDS(paste0(out.dir,"pre-transplant_post-transplant_integration_replicate2.rds"))


# Fig2I-M - Putative TF regulators for pre- and post-transplant clusters -------
#Notes: 
# TF candidate analysis where comparison is only on the pre- or post-transplant
# groups.

# Select clustering resolution
resolution <- "integrated_snn_res.0.4"

# Set identity class to resolution 0.4 annotations
Idents(object = ilc.integrated) <- resolution

# Store cluster names to iterate through
clusters.names <- levels(ilc.integrated)

# Store number of motifs at gained regions of chromatin accessibility
num.motifs.at.gained <- list()

# Store number of putative regulators identified for each cluster 
num.putative.regulators <- list()

#
cluster.labels <- c("Pre2", "Pre3", "Pre1",
                    "Post1", "Post2", "Post3")

colors <- c("#63C7C5", "#A77F48", "#005BA3",
            "#238755", "#88B571", "#57813D")


for (i in 1:length(clusters.names)) {
  cluster = clusters.names[i]
  
  ## Putative regulators of pre-transplant clusters ----------------------------
  if (cluster == 0 | cluster == 1 | cluster == 2) {
    
    # Differential RNA analysis 
    markers_rna <- presto:::wilcoxauc.Seurat(X = ilc.integrated,
                                             group_by = resolution,
                                             assay = 'data',
                                             seurat_assay = 'SCT',
                                             groups_use = c(cluster, "3","4","5"))
    
    # Differential motif analysis 
    markers_motifs <- presto:::wilcoxauc.Seurat(X = ilc.integrated,
                                                group_by = resolution,
                                                assay = 'data',
                                                seurat_assay = 'chromvar',
                                                groups_use = c(cluster, "3","4","5"))
    
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
    # https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis.html
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
    
    output.filename = paste0(out.dir,"Fig2l_", cluster.labels[i], "_cluster", cluster, "_v_postClusters_candidate_TF.pdf")
    
    ggsave(file = output.filename,height = 4.6, width = 3, device = "pdf")
    
  }else if (cluster == 3 | cluster == 4 | cluster == 5) {
    ## Putative regulators of post-transplant clusters ---------------------------  
    # Differential RNA analysis 
    markers_rna <- presto:::wilcoxauc.Seurat(X = ilc.integrated,
                                             group_by = resolution,
                                             assay = 'data',
                                             seurat_assay = 'SCT',
                                             groups_use = c(cluster, "0","1","2"))
    
    # Differential motif analysis 
    markers_motifs <- presto:::wilcoxauc.Seurat(X = ilc.integrated,
                                                group_by = resolution,
                                                assay = 'data',
                                                seurat_assay = 'chromvar',
                                                groups_use = c(cluster, "0","1","2"))
    
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
    # Source of function:
    # https://satijalab.org/seurat/articles/weighted_nearest_neighbor_analysis.html
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
    
    
    output.filename = paste0(out.dir,"Fig2m_", cluster.labels[i], "_cluster", cluster,"_v_preClusters_candidate_TF.pdf")
    
    ggsave(file = output.filename,height = 4.6, width = 3, device = "pdf")
    
    ## Volcano plot of differential motifs ---------------------------------------
    
    # Direction of Log2Fold change is relative to a particular cluster 
    # Here I've selected either the pre- or post- cluster that is be tested
    # against the entirety of post- or pre- samples respectively
    motif.data <- markers_motifs %>% 
      dplyr::filter(motif.group == cluster) %>% 
      arrange(desc(motif.logFC))
    
    # Add column that specifies the color of significant or non-significant 
    # motifs
    motif.volcano.data <- motif.data %>% 
      mutate(sig.color  = ifelse(-log10(motif.padj) > -log10(0.05), '#1F9FE1', '#E0601E')) %>%
      arrange(sig.color)
    
    # Volcano Plot - highlighting specific ILC1 and ILC2 genes
    ggplot(motif.volcano.data, aes(x=motif.logFC, y=-log10(motif.padj),  col = sig.color, label = gene)) + 
      geom_point() +
      scale_color_identity() +
      geom_text_repel(box.padding = .5,
                      size=3,
                      segment.size = 0.25,
                      force_pull=1,
                      max.overlaps=50,
                      min.segment.length= 0.01,
                      data=subset(motif.volcano.data, (gene == "GATA3" | gene == "TBX21" | gene == "BATF" |
                                                         gene == "MAF" | gene == "FOSB" | gene == "BACH1" |
                                                         gene == "EOMES" | gene == "RBPJL" | gene == "FOXO3" |
                                                         gene == "TCF3" | gene == "SMAD2" | gene == "BATF" |
                                                         gene == "NFATC4" | gene == "NR4A1" | gene == "CREB1" |
                                                         gene == "MECOM" | gene == "IRF8" | gene == "NR4A2" |
                                                         gene == "POU2F2" | gene == "SNAI3")),
                      colour = 'black') +
      theme_classic() + 
      theme(legend.position="none") +
      ylab("-log10(p-value)") +
      xlab(paste0("Motif Enrichment \n log2(",cluster.labels[i], "/ Pre 1-3)"))
    
    # File name parameters
    if (cluster == 3) {
      fig = "Fig2i_"
    }else if( cluster == 4){
      fig = "Fig2j_"
    }else{
      fig = "Fig2k_"
    }
    
    output.filename.vol = paste0(out.dir, fig, cluster.labels[i],"_post_motif_volcano_cluster_",cluster, ".pdf")
    
    # Save Volcano plot
    ggsave(file = output.filename.vol,
           height = 4.5, width = 4.5, device = "pdf")
  }
  
  
  ## Number of differential motifs at gained regions ---------------------------
  
  # Number of differential motifs at gained regions of chromatin accessibility
  num.motifs.cluster <- markers_motifs %>% 
    dplyr::filter(motif.group == cluster & motif.logFC > 0 & -log10(motif.padj) > -log10(0.05)) %>%
    arrange(desc(motif.logFC)) %>% nrow()
  
  num.motifs.at.gained <- append(num.motifs.at.gained, list(num.motifs.cluster))
  
}


# Save Session Info ------------------------------------------------------------

writeLines(capture.output(sessionInfo()), "./sessionInfo/Figure_2i_m_sessionInfo.txt")
