# Performs differential chromatin accessibility analysis across the 6 clusters. 
# Each of the pre clusters are compared to all of the post-transplant clusters, 
#and each of the post-transplant clusters are compared to all of the 
#pre-transplantv clusters.
#
# Generates plot for Figure 2h


library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)

# Output Directory -------------------------------------------------------------

out.dir <- "./Figure_2/"


# Read Integrated Seurat Object ------------------------------------------------

ilc.integrated <- readRDS(paste0(out.dir,"pre-transplant_post-transplant_integration_replicate2.rds"))


# FIG. 2h ----------------------------------------------------------------------

# Note: Variable names have been carried over from the 
# RNA abundance analysis (Figure 2h)

## Differential Chromatin Accessibility analysis --------------------------------

# Select clusters for a given resolution value
Idents(object = ilc.integrated) <- "integrated_snn_res.0.4"

# Store differential regions of chromatin accessibility for each cluster
diff_gene_list <- list()

# pre-transplant differential ATAC peaks relative to post-transplant cells
pre.transplant.clusters <- c(1,0,2)

for (pre.clust in pre.transplant.clusters) {
  
  # Identify differential ATAC peaks
  diff.genes <- FindMarkers(ilc.integrated,
                            ident.1 = pre.clust,
                            ident.2 = c(3,4,5),
                            min.pct = 0.25,
                            assay = "ATAC",
                            slot = "data")  
  
  # Select significant ATAC peaks with FC > 0
  top.clusterX.genes <- diff.genes %>%
    dplyr::filter(p_val_adj < 0.05 & avg_log2FC > 0) %>%
    arrange(desc(avg_log2FC))
  
  diff_gene_list <- append(diff_gene_list,list(row.names(top.clusterX.genes)))

} 


# post-transplant differential AAC peaks relative to all pre-transplant cells
post.transplant.clusters <- c(3,4,5)

for (post.clust in post.transplant.clusters) {
  
  # Identify differential genes
  diff.genes <- FindMarkers(ilc.integrated,
                            ident.1 = post.clust,
                            ident.2 = c(0,1,2),
                            min.pct = 0.25,
                            assay = "ATAC",
                            slot = "data")  
  
  # Select significant ATAC peaks with FC > 0
  top.clusterX.genes <- diff.genes %>%
    dplyr::filter(p_val_adj < 0.05 & avg_log2FC > 0) %>%
    arrange(desc(avg_log2FC))
  
  diff_gene_list <- append(diff_gene_list, list(row.names(top.clusterX.genes)))


} 

## Generate Union Set of differential genes ------------------------------------

# Clusters to cycle through
cluster.names <- levels(ilc.integrated)

#Need vector of cell IDs for each cluster 
cluster.cellIDs <- list()

#List of matricies of average matricies 
avg_gene_exp <- list()

# List of avg gene dataframes
gene_set_avg <- list()

for (gene_set in diff_gene_list) {
  
  for (cluster in cluster.names) {
    
    # Isolate the cellIDs for each cluster
    cell_IDs <- row.names(ilc.integrated@meta.data %>% dplyr::filter(integrated_snn_res.0.4 == cluster))
    
    # Normalized ATAC data
    rna_data <- ilc.integrated@assays[["ATAC"]]@data
    
    # Sub-matrix of normalized RNA signal
    sub_rna_data <- rna_data[gene_set, cell_IDs]
    
    # Vector of avg ATAC signal
    sub_rna_avg <- rowMeans(sub_rna_data)
    
    # Store avg ATAC signal (across cells) for each cluster
    avg_gene_exp <- append(avg_gene_exp, list(sub_rna_avg))
    
  }
  
  #Concat average ATAC signal vectors into data frame for each ATAC peak set
  gene_set_avg_df <-as.data.frame(do.call(cbind, avg_gene_exp))
  
  # Save avg. ATAC signal for each ATAC peak set for each cluster
  gene_set_avg <- append(gene_set_avg, list(gene_set_avg_df))
  
  # Reset the list holding the average score for each gene set
  avg_gene_exp <- list()
  
}

# Full Matrix of ATAC peaks 
diff_gene_clust_data <- do.call(rbind, gene_set_avg)

## Complex Heatmap -------------------------------------------------------------

# Mean Center the rows 
df <- t(scale(t(data.matrix(diff_gene_clust_data)), scale = FALSE))

# Remove mean attribute from matrix - Not necessary 
attr(df, "scaled:center") <- NULL

cluster_names = c("Pre1","Pre2","Pre3",
                  "Post1","Post2","Post3")

row_split = rep(cluster_names, times = lengths(diff_gene_list))

#order row split with factor
row_split <- factor(row_split, levels=cluster_names)

# Add column names
colnames(df) <- cluster_names

# Colors for Heatmap
atac_color = colorRamp2(c(-0.4, 0, 0.4), c("#175e7a", "white", "#70583a"))
atac_color(seq(-0.4, 0.4))

pdf(paste0(out.dir, "Fig2h_pre_post_cluster_ATAC_heatmap.pdf"),
    width=3,
    height=3.5)

Heatmap(df,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        row_names_gp = gpar(fontsize = 0),
        name = "Norm.\nGene Exp.",
        left_annotation = HeatmapAnnotation(ATAC_Peak_Sets = row_split,
                                            col = list(ATAC_Peak_Sets = c("Pre1" = "#005ba3",
                                                                    "Pre2" = "#a77f48",
                                                                    "Pre3" = "#63c7c5",
                                                                    "Post1" = "#238755",
                                                                    "Post2" = "#88b571",
                                                                    "Post3" = "#57813d")),
                                            which = 'row'),
        row_split = row_split,
        row_title = NULL,
        col = atac_color)


dev.off()

# Save Session Info ------------------------------------------------------------

writeLines(capture.output(sessionInfo()), "./sessionInfo/Figure_2h_sessionInfo.txt")
