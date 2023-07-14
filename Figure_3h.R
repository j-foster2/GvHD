# Generates Figure 3h


library(tidyverse)
library(Seurat)
library(Signac)
library(stringr)
library(ggplot2)


# Output Directory -------------------------------------------------------------
out.dir <- "./Figure_3/"


# Read Replicate 1 Data - Seurat Object ----------------------------------------
ilc.integrated <- readRDS(file = "./Figure_3/inVitro_ILC2_pcILC2_integrated_rep1.rds")

# Add gene activity ------------------------------------------------------------

# quantify gene activity
gene.activities <- GeneActivity(ilc.integrated)

# add gene activities as a new assay
ilc.integrated[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)

# normalize gene activities - per multimodal reference mapping
DefaultAssay(ilc.integrated) <- "ACTIVITY"

ilc.integrated <- SCTransform(ilc.integrated, assay = "ACTIVITY", new.assay.name = "ACTIVITY_sct", verbose = FALSE)


# ILC2 associated genes - GBA --------------------------------------------------
ilc2.genes <- read.table("./referenceData/ILC2_GBA_genes.bed")

# Bulk GBA
ilc2.genes <- read.table("./referenceData/GBA_ILC2_genes_bulk.txt",
                         skip = 0)

colnames(ilc2.genes) <- "geneNames"

# remove ILC2 genes not in the activity matrix
ilc2.genes <- rownames(ilc.integrated@assays[["ACTIVITY_sct"]]@data)[rownames(ilc.integrated@assays[["ACTIVITY_sct"]]@data) %in% ilc2.genes$geneNames]

ilc2.gba.sig <- colMeans(ilc.integrated@assays[["ACTIVITY_sct"]]@data[ilc2.genes,])

ilc.integrated@meta.data['ilc2.gba.sig'] <- ilc2.gba.sig


# Update the orig ident 
ilc.integrated@meta.data['clust_ilc_classification'] <- with(ilc.integrated@meta.data, ifelse(seurat_clusters %in% c(2,3,4,5,6,7), 'ILC2',
                                                                                              ifelse(seurat_clusters %in% c(0,1,8), 'pcILC2', 'ILC1-like')))


# Plot the Distribution of GBA associated genes --------------------------------
ilc.integrated$clust_ilc_classification <- factor(ilc.integrated$clust_ilc_classification, levels = c("ILC2","pcILC2",'ILC1-like'))

# ILC2 Plot 
fig3h <- VlnPlot(object = ilc.integrated,
        features = c("ilc2.gba.sig"),
        group.by = "clust_ilc_classification", 
        pt.size = 0,
        same.y.lims = T,
        cols = c("#2C5D73","#F0972B","#6A4C87")) +
  geom_boxplot(width=0.1)  +
  stat_summary(fun=mean, geom="point", size=2, color="red") + 
  theme(legend.position="none") +
  ggtitle("")

ggsave(fig3h,
       file = paste0(out.dir,"Fig3h_ILC2_pcILC2_ILC1-like_avgChromAccess_ILC2_GBA_genes.pdf"),
       height = 5.5,
       width = 4.5,
       device = "pdf")


# Wilcoxon test ----------------------------------------------------------------

ilc2.chrom.signal.at.ilc2 <- ilc.integrated@meta.data %>% dplyr::filter(clust_ilc_classification == "ILC2") %>% dplyr::select(ilc2.gba.sig)

pc.ilc2.chrom.signal.at.ilc2 <- ilc.integrated@meta.data %>% dplyr::filter(clust_ilc_classification == "pcILC2") %>% dplyr::select(ilc2.gba.sig)

ilc1.like.chrom.signal.at.ilc2 <- ilc.integrated@meta.data %>% dplyr::filter(clust_ilc_classification == "ILC1-like") %>% dplyr::select(ilc2.gba.sig)

# ILC2 v. pcILC2
ilc2_pcilc2_wilcox <- wilcox.test(x = ilc2.chrom.signal.at.ilc2$ilc2.gba.sig,
            y = pc.ilc2.chrom.signal.at.ilc2$ilc2.gba.sig,
            alternative = "greater")

# ILC2 v. ILC1-like
ilc2_ilc1_like_wilcox <- wilcox.test(x = ilc2.chrom.signal.at.ilc2$ilc2.gba.sig,
            y = ilc1.like.chrom.signal.at.ilc2$ilc2.gba.sig,
            alternative = "greater")


# Save Mann Whitney U Test Results for Fig3 ------------------------------------

write("Mann Whitney U Test Results for Fig3h:\n",
      file= paste0(out.dir, "Fig3h_stats.txt"))

write("ILC2 v. pcILC2:",
      file= paste0(out.dir, "Fig3h_stats.txt"),
      append=TRUE)

write(capture.output(ilc2_pcilc2_wilcox),
      file= paste0(out.dir, "Fig3h_stats.txt"),
      append=TRUE)

write("ILC2 v. ILC1-like:",
      file= paste0(out.dir, "Fig3h_stats.txt"),
      append=TRUE)

write(capture.output(ilc2_ilc1_like_wilcox),
      file= paste0(out.dir, "Fig3h_stats.txt"),
      append=TRUE)


# Save Session Info ------------------------------------------------------------
writeLines(capture.output(sessionInfo()), "./sessionInfo/Figure_3h_sessionInfo.txt")

