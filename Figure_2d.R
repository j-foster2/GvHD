# Performs differential expression analysis across the 6 clusters. Each of the
# pre clusters are compared to all of the post-transplant clusters, and each 
# of the post-transplant clusters are compared to all of the pre-transplant 
# clusters.
#
# Generates plot for Figure 2d


library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)

# Output Directory -------------------------------------------------------------

out.dir <- "./Figure_2/"


# Read Integrated Seurat Object ------------------------------------------------

ilc.integrated <- readRDS(paste0(out.dir,"pre-transplant_post-transplant_integration_replicate2.rds"))


# FIG. 2d - Heatmap - Average Gene Expression ----------------------------------

# Differential gene expression analysis of clusters 

Idents(object = ilc.integrated) <- "integrated_snn_res.0.4"

# Store differential genes for each cluster
diff_gene_list <- list()

# pre-transplant differential genes relative to post-transplant cells
pre.transplant.clusters <- c(2,0,1)

for (pre.clust in pre.transplant.clusters) {
  
  # Identify differential genes
  diff.genes <- FindMarkers(ilc.integrated,
                            ident.1 = pre.clust,
                            ident.2 = c(3,4,5),
                            min.pct = 0.25,
                            assay = "SCT")  
  
  # Select significant genes with FC > 0
  top.clusterX.genes <- diff.genes %>%
    dplyr::filter(p_val_adj < 0.05 & avg_log2FC > 0) %>%
    arrange(desc(avg_log2FC))
  
  diff_gene_list <- append(diff_gene_list,list(row.names(top.clusterX.genes)))
  
  # Write text file of differential genes for each cluster
  write.table(top.clusterX.genes, 
              paste0(out.dir, "ILC2_transplant_differential_genes_cluster_",pre.clust,".txt"),
              sep = '\t',
              row.names = T,
              quote = F)
  
}

# post-transplant differential genes relative to all pre-transplant cells
post.transplant.clusters <- c(3,4,5)

for (post.clust in post.transplant.clusters) {
  
  # Identify differential genes
  diff.genes <- FindMarkers(ilc.integrated,
                            ident.1 = post.clust,
                            ident.2 = c(0,1,2),
                            min.pct = 0.25,
                            assay = "SCT")  
  
  # Select significant genes with FC > 0
  top.clusterX.genes <- diff.genes %>%
    dplyr::filter(p_val_adj < 0.05 & avg_log2FC > 0) %>%
    arrange(desc(avg_log2FC))
  
  diff_gene_list <- append(diff_gene_list, list(row.names(top.clusterX.genes)))
  
  # Write text file of differential genes for each cluster
  write.table(top.clusterX.genes,paste0(out.dir, "ILC2_transplant_differential_genes_cluster_",post.clust,".txt"),
              sep = '\t',
              row.names = T,
              quote = F)
}

# Store list of differential genes for each cluster
saveRDS(diff_gene_list, file = paste0(out.dir, "Fig2d_gene_sets.rds"))

# Generate Union Set of differential genes

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
    
    # Normalized RNA data
    rna_data <- ilc.integrated@assays[["SCT"]]@data
    
    # Submatrix of normalized RNA signal
    sub_rna_data <- rna_data[gene_set, cell_IDs]
    
    # Vector of avg gene expression
    sub_rna_avg <- rowMeans(sub_rna_data)
    
    # Store avg gene values (across cells) for each cluster
    avg_gene_exp <- append(avg_gene_exp, list(sub_rna_avg))
    
  }
  
  #Concat avgerage gene vectors into dataframe for each gene set
  gene_set_avg_df <-as.data.frame(do.call(cbind, avg_gene_exp))
  
  # Save the 
  gene_set_avg <- append(gene_set_avg, list(gene_set_avg_df))
  
  # Reset the list holding the average score for each gene set
  avg_gene_exp <- list()
  
}

# Highlight pre-transplant shared genes in across pre groups
pre_shared_genes <- intersect(intersect(rownames(gene_set_avg[[1]]),rownames(gene_set_avg[[2]])),rownames(gene_set_avg[[3]]))

# Highlight post-transplant shared genes in across post groups
post_shared_genes <- intersect(intersect(rownames(gene_set_avg[[4]]),rownames(gene_set_avg[[5]])),rownames(gene_set_avg[[6]]))

# Copy gene names to column (needed to identified duplicates)
gene_set_avg_keepNames <- list()
for (df in gene_set_avg) {
  
  df$geneName <- rownames(df)
  
  gene_set_avg_keepNames <- append(gene_set_avg_keepNames, list(df))
  
}

# Full Matrix of genes 
diff_gene_clust_data <- do.call(rbind, gene_set_avg_keepNames)

# pre-transplant shared genes boolean of shared or not
pre_shared_genes <- c("pre_not_shared","pre_shared")[1 + (diff_gene_clust_data[,'geneName'] %in% pre_shared_genes)]

# post-transplant shared genes boolean of shared or not
post_shared_genes <- c("post_not_shared","post_shared")[1 + (diff_gene_clust_data[,'geneName'] %in% post_shared_genes)]

# Drop geneName column 
diff_gene_clust_data <- diff_gene_clust_data[-c(7)]

# Complex Heatmap

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

pdf(paste0(out.dir, "Fig2d_pre_post_cluster_heatmap.pdf"),
    width=4,
    height=4.25)

Heatmap(df,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        row_names_gp = gpar(fontsize = 0),
        name = "Norm.\nGene Exp.",
        left_annotation = HeatmapAnnotation(gene_sets = row_split,
                                            col = list(gene_sets = c("Pre1" = "#63c7c5",
                                                                    "Pre2" = "#a77f48",
                                                                    "Pre3" = "#005ba3",
                                                                    "Post1" = "#238755",
                                                                    "Post2" = "#88b571",
                                                                    "Post3" = "#57813d")),
                                            which = 'row'),
        right_annotation = HeatmapAnnotation(Pre_Intersect = pre_shared_genes,
                                             Post_Intersect = post_shared_genes,
                                             col = list(Pre_Intersect = c("pre_shared" = "#000000",
                                                                          "pre_not_shared" = "#FFFFFF"),
                                                        Post_Intersect = c("post_shared" = "#000000",
                                                                           "post_not_shared" = "#FFFFFF")),
                                             which = 'row',
                                             show_legend = FALSE),
        
        row_split = row_split,
        row_title = NULL)


dev.off()

# Save Session Info ------------------------------------------------------------

writeLines(capture.output(sessionInfo()), "./sessionInfo/Figure_2d_sessionInfo.txt")
