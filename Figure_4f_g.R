# Generates Figures 4f-g


library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(cowplot)
library(data.table)
library(tidyverse)
library(GenomicRanges)
library(biomaRt)
library(DESeq2)
library(gprofiler2)


# Output Directory -------------------------------------------------------------
out.dir <- "./Figure_4/"


# Read merged human scATAC data ------------------------------------------------
combined <- readRDS(file = paste0(out.dir,"human_scATAC_Merged_plusPseudoExp_SeuratObject.rds"))


# Create Count matrix for each sample ------------------------------------------
sum_expression_data <- AggregateExpression(object = combined,
                                           assays = "ACTIVITY",
                                           slot = "counts",
                                           group.by = "dataset")

sum_expression_data_df <- as.data.frame(sum_expression_data)

# Update column names
sampleIDs <- c("ILC2","Healthy1","Healthy2",
               "PreNoGvHD1", "PreNoGvHD2", "PreNoGvHD3",
               "PreGvHD1", "PreGvHD2", "PreGvHD3",
               "PostNoGvHD1", "PostNoGvHD2", "PostNoGvHD3",
               "PostGvHD1", "PostGvHD2", "PostGvHD3")

colnames(sum_expression_data_df) <- sampleIDs

# Select pre and post-transplant samples
sum_expression_data_prePost_df <- sum_expression_data_df[-c(1:3)]


# DESeq2 analysis --------------------------------------------------------------

# Create coldata file
coldata_df <- as.data.frame(c("PreNoGvHD", "PreNoGvHD", "PreNoGvHD",
                              "PreGvHD", "PreGvHD", "PreGvHD",
                              "PostNoGvHD", "PostNoGvHD", "PostNoGvHD",
                              "PostGvHD", "PostGvHD", "PostGvHD"))

rownames(coldata_df) <- c("PreNoGvHD1", "PreNoGvHD2", "PreNoGvHD3",
                          "PreGvHD1", "PreGvHD2", "PreGvHD3",
                          "PostNoGvHD1", "PostNoGvHD2", "PostNoGvHD3",
                          "PostGvHD1", "PostGvHD2", "PostGvHD3")

coldata_df$samples <- c("C","D","E",
                        "F","G", "H",
                        "K","L","M",
                        "N","O","P")

colnames(coldata_df) <- c("codition", "sample")

coldata_df$pre_post <- c(
  "pre","pre","pre",
  "pre","pre", "pre",
  "post","post","post",
  "post","post","post")

# Set Pre- Post- column as factor
coldata_df$pre_post <-factor(coldata_df$pre_post, levels = c("post","pre"))

coldata_df$codition <- factor(coldata_df$codition)

coldata_df$sample <- factor(coldata_df$sample)

dds <- DESeqDataSetFromMatrix(sum_expression_data_prePost_df, coldata_df, ~pre_post)

dds = DESeq(dds)

deseq2_results = results(dds)

results_df <- as.data.frame(deseq2_results)

# Log2 transform the base Mean 
results_df$log2baseMean <- log2(results_df$baseMean)

# Fig. 4f - Volcano Plot -------------------------------------------------------

results_df$neg.log10padj <- -log10(results_df$padj)
pseudo_volcano_plot <-results_df %>% dplyr::arrange(desc(padj)) %>%
  ggplot(aes(x = log2FoldChange, y = neg.log10padj, colour = padj < 0.05)) +
  geom_point() +
  theme_classic() +
  scale_colour_manual(values = setNames(c('#005a9e','#7f7f7f'),c(T, F))) +
  theme(legend.position="none") +
  aes(stroke = 0.25) +
  xlab("ATAC signal gene bodies \n log2(Post/Pre-transplant)") +
  ylab("-log10(padj)")

ggsave(pseudo_volcano_plot,
       file = paste0(out.dir, "Fig4f_human_scATAC_pre_post_pseudogene_diff_volcanoPlot.pdf"),
       height = 4.4,
       width = 5.5 ,
       device = "pdf")


# Pre-transplant differential genes --------------------------------------------
baseMean_val = 1000

dif_genes_pre <- results_df %>%
  dplyr::filter(baseMean > baseMean_val & log2FoldChange < 0, padj < 0.05) %>%
  dplyr::arrange(log2FoldChange)

write.table(rownames(dif_genes_pre),
            file = paste0(out.dir, "human_pre_transplant_diff_genes.txt"),
            sep = '\t',
            col.names = F,
            quote = FALSE,
            row.names = FALSE)

# Post-transplant differential genes -------------------------------------------
dif_genes_post <- results_df %>% dplyr::filter(baseMean > baseMean_val & log2FoldChange > 0, padj < 0.05) %>% dplyr::arrange(desc(log2FoldChange))

write.table(rownames(dif_genes_post),
            file = paste0(out.dir, "human_post_transplant_diff_genes.txt"),
            sep = '\t',
            col.names = F,
            quote = FALSE,
            row.names = FALSE)

# Pre-transplant g:Profiler analysis -------------------------------------------

# List of differential genes
pre_diff_genes <- row.names(dif_genes_pre)

pre_gost <- gost(query = pre_diff_genes, 
                 organism = "hsapiens", ordered_query = TRUE, 
                 multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                 measure_underrepresentation = FALSE, evcodes = FALSE, 
                 user_threshold = 0.05, correction_method = "g_SCS", 
                 domain_scope = "annotated", custom_bg = NULL, 
                 numeric_ns = "", sources = NULL, as_short_link = FALSE)

pre_gost$result$neglog10padj <- -log10(pre_gost$result$p_value)

pre_goterms_plot <- head(pre_gost$result, n =5) %>%
  ggplot(aes(x = neglog10padj, y = reorder(term_name, neglog10padj))) +
  geom_bar(stat = "identity", width = 0.75) +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) + 
  ylab("GO:BP Term") + 
  theme_classic() 

ggsave(pre_goterms_plot,
       file = paste0(out.dir, "Fig4g_human_scATAC_pre_diff_GoTerms.pdf"),
       height = 4.4,
       width = 5.5 ,
       device = "pdf")

# Post-transplant g:Profiler analysis ------------------------------------------

# List of differential genes
post_diff_genes <- row.names(dif_genes_post)

post_gost <- gost(query = post_diff_genes, 
                  organism = "hsapiens", ordered_query = TRUE, 
                  multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                  measure_underrepresentation = FALSE, evcodes = FALSE, 
                  user_threshold = 0.05, correction_method = "g_SCS", 
                  domain_scope = "annotated", custom_bg = NULL, 
                  numeric_ns = "", sources = NULL, as_short_link = FALSE)

post_gost$result$neglog10padj <- -log10(post_gost$result$p_value)

post_goterms_plot <- head(post_gost$result, n =5) %>%
  ggplot(aes(x = neglog10padj, y = reorder(term_name, neglog10padj))) +
  geom_bar(stat = "identity", width = 0.75) +
  scale_y_discrete(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) + 
  ylab("GO:BP Term") + 
  theme_classic() 

ggsave(post_goterms_plot,
       file = paste0(out.dir, "Fig4g_human_scATAC_post_diff_GoTerms.pdf"),
       height = 4.4,
       width = 5.5 ,
       device = "pdf")

# Save Session Info ------------------------------------------------------------
writeLines(capture.output(sessionInfo()), "./sessionInfo/Figure_4f_g_sessionInfo.txt")
