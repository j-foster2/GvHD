# Calculate ATAC signal over gene bodies

library(tidyverse)
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(GenomicRanges)


# Output Directory -------------------------------------------------------------
out.dir <- "./Figure_4/"


# Read combined human scATAC data ----------------------------------------------
combined <- readRDS(paste0(out.dir,"human_scATAC_Merged_SeuratObject.rds"))


# Preliminary Analysis ---------------------------------------------------------
combined <- RunTFIDF(combined)

combined <- FindTopFeatures(combined, min.cutoff = 20)

combined <- RunSVD(combined)

combined <- RunUMAP(combined, dims = 2:50, reduction = 'lsi')


# Annotating ChromatinAssayObject with genome information ----------------------
# Reference GRanges for hg38
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

combined@assays[["ATAC"]]@annotation <- annotations
genome(combined) <- "hg38"


# Identifying anchors between scRNA-seq and scATAC-seq datasets ----------------

# quantify gene activity
gene.activities <- GeneActivity(combined)

# add gene activities as a new assay
combined[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)

# normalize gene activities - per multimodal reference mapping
DefaultAssay(combined) <- "ACTIVITY"

combined <- SCTransform(combined,assay = "ACTIVITY", verbose = FALSE)


# Save Merged + PsuedoExpression Seurat Object ---------------------------------
saveRDS(combined, file = paste0(out.dir,"human_scATAC_Merged_plusPseudoExp_SeuratObject.rds"))


# Save Session Info ------------------------------------------------------------
writeLines(capture.output(sessionInfo()), "./sessionInfo/addPseudoExpression_integrated_human_scATAC_sessionInfo.txt")
