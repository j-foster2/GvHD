# Integration of the second replicate of single cell Multiome data 
#
# Used for Figure 2b, c

# Memory Requirement: 

library(JASPAR2020)
library(TFBSTools)
library(Seurat)
library(EnsDb.Mmusculus.v79)
library(GenomicAlignments)
library(Signac)
library(patchwork)
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggrepel)
library(gprofiler2)
library(forcats)
library(clustree)
library(gplots)
library(UpSetR)
library(chromVAR)
library(motifmatchr)
library(BSgenome.Mmusculus.UCSC.mm10)
library(circlize)
library(ComplexHeatmap)

# setwd("/proj/dllab/jfoster/serody_project/GitHub/GvHD/")

# Output Directory -------------------------------------------------------------

out.dir <- "./Figure_2/"


# Load Motif Database ----------------------------------------------------------

pwm_set <- getMatrixSet(x = JASPAR2020,
                        opts = list(collection = "CORE",
                                    tax_group = 'vertebrates',
                                    all_versions = FALSE))


# Reference GRanges for mm10 ---------------------------------------------------

annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"


# Load Multiome data - HDF5 format ---------------------------------------------

#Load Pre-transplant and post-transplant data 
pre_transplant.data <- Read10X_h5("./data/ilc2_pre-transplant_rep2_filtered_feature_bc_matrix.h5")

post_transplant.data <- Read10X_h5("./data/ILC2_post-transplant_rep2_filtered_feature_bc_matrix.h5")


# Read Pre-transplant ILC2 data - (including peak identification) --------------

# Extract RNA and ATAC data
pre_transplant.rna_counts = pre_transplant.data$`Gene Expression`

pre_transplant.atac_counts <- pre_transplant.data$Peaks

# Create Seurat object for pre-transplant data
pre_transplant = CreateSeuratObject(counts = pre_transplant.rna_counts, project = "pre_transplant")

# Add ATAC data to Seurat object
grange.counts <- StringToGRanges(rownames(pre_transplant.atac_counts), sep = c(":", "-"))

grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)

pre_transplant.atac_counts <- pre_transplant.atac_counts[as.vector(grange.use), ]

# Store snATAC Fragment file
pre_transplant_frag.file <- "./data/ilc2_pre-transplant_rep2_atac_fragments.tsv.gz"

chrom_assay <- CreateChromatinAssay(
  counts = pre_transplant.atac_counts,
  sep = c(":", "-"),
  genome = 'mm10',
  fragments = pre_transplant_frag.file,
  min.cells = 10,
  annotation = annotations
)

pre_transplant[["ATAC"]] <- chrom_assay


# Read Post-transplant ILC2 data - (including peak identification) -------------

# Extract RNA and ATAC data
post_transplant.rna_counts = post_transplant.data$`Gene Expression`

post_transplant.atac_counts <- post_transplant.data$Peaks

# Create Seurat object for GFP data
post_transplant = CreateSeuratObject(counts = post_transplant.rna_counts,
                                     project = "post_transplant")

# Add ATAC data to Seurat Object
grange.counts <- StringToGRanges(rownames(post_transplant.atac_counts), sep = c(":", "-"))

grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)

post_transplant.atac_counts <- post_transplant.atac_counts[as.vector(grange.use), ]


post_transplant.frag.file <- "./data/ILC2_post-transplant_rep2_atac_fragments.tsv.gz"

chrom_assay <- CreateChromatinAssay(
  counts = post_transplant.atac_counts,
  sep = c(":", "-"),
  genome = 'mm10',
  fragments = post_transplant.frag.file,
  min.cells = 10,
  annotation = annotations
)

post_transplant[["ATAC"]] <- chrom_assay


# QC Filtering -----------------------------------------------------------------

pre_transplant[["percent.mt"]] <- PercentageFeatureSet(pre_transplant, pattern = "^mt-")

post_transplant[["percent.mt"]] <- PercentageFeatureSet(post_transplant, pattern = "^mt-")

# Filter pre-transplant cells
pre_transplant <- subset(
  x = pre_transplant,
  subset = nCount_ATAC < 7e4 &
    nCount_ATAC > 5e3 &
    nCount_RNA < 25000 &
    nCount_RNA > 1000 &
    percent.mt < 20
)

# Filter post-transplant cells
post_transplant <- subset(
  x = post_transplant,
  subset = nCount_ATAC < 7e4 &
    nCount_ATAC > 5e3 &
    nCount_RNA < 25000 &
    nCount_RNA > 1000 &
    percent.mt < 20
)


# Integrate pre- and post-transplant data --------------------------------------

DefaultAssay(pre_transplant) <- "RNA"

DefaultAssay(post_transplant) <- "RNA"

# List of Seurat Objects
ILC.list = list(pre_transplant, post_transplant)

# Normalization and identification of variable features
ILC.list = lapply(X = ILC.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# Identify features that are variable across data sets
features = SelectIntegrationFeatures(object.list = ILC.list)

immune.anchors = FindIntegrationAnchors(object.list = ILC.list,
                                        anchor.features = features)

# Integrate data using RNA data
ilc.integrated = IntegrateData(anchorset = immune.anchors)

# Run Seurat workflow for Integrated data
DefaultAssay(ilc.integrated) <- "integrated"

# TODO scale data - describe mroe
ilc.integrated = ScaleData(ilc.integrated, verbose = FALSE)

# Dimension reduction with PCA
ilc.integrated = RunPCA(ilc.integrated, npcs = 30, verbose = FALSE)

# Dimension reduction and visualization of data with UMAP
ilc.integrated = RunUMAP(ilc.integrated,
                         reduction = "pca",
                         dims = 1:30)

# Compute nearest neighbors
ilc.integrated = FindNeighbors(ilc.integrated, reduction = "pca", dims = 1:30)

# Cluster identification at several resolutions
ilc.integrated = FindClusters(ilc.integrated, resolution = c(0, 0.05, 0.1, 0.2,
                                                             0.3, 0.4, 0.5,
                                                             0.6, 0.7, 0.8,
                                                             0.9, 1.0))


# Evaluate impact of altering resolution for clustering ------------------------

ilc.clusttree <- clustree(ilc.integrated,
                          prefix = "integrated_snn_res.")

# Store the clusttree that shows there are a handful of cells that
# are always in a separate cluster. This holds true even at low
# resolution values.
ilc.integrated@misc <- append(ilc.integrated@misc, list(ilc.clusttree))

names(ilc.integrated@misc) <- c("clusterOptimizationPlot")

#View plot
ilc.integrated@misc$clusterOptimizationPlot


# Drop two cells that cluster on their own at multiple resolution values -------

cells.to.remove <- c(row.names(ilc.integrated@meta.data %>%
                                 dplyr::filter(ilc.integrated@meta.data['integrated_snn_res.0.05'] == 2)))

ilc.integrated <- ilc.integrated[,!colnames(ilc.integrated) %in% cells.to.remove]


# RNA analysis (Seurat workflow) -----------------------------------------------

DefaultAssay(ilc.integrated) <- "RNA"

ilc.integrated <- SCTransform(ilc.integrated, verbose = FALSE) %>%
  RunPCA() %>%
  RunUMAP(dims = 1:50,
          reduction.name = 'umap.rna',
          reduction.key = 'rnaUMAP_')


# ATAC analysis (Seurat workflow) ----------------------------------------------

DefaultAssay(ilc.integrated) <- "ATAC"

# Term-frequency inverse-document-frequency Normalization
ilc.integrated <- RunTFIDF(ilc.integrated)

ilc.integrated <- FindTopFeatures(ilc.integrated, min.cutoff = 'q0')

ilc.integrated <- RunSVD(ilc.integrated)

ilc.integrated <- RunUMAP(ilc.integrated,
                          reduction = 'lsi',
                          dims = 2:50,
                          reduction.name = "umap.atac",
                          reduction.key = "atacUMAP_")

# Identification of motifs
DefaultAssay(ilc.integrated) <- "ATAC"

motif.matrix <- CreateMotifMatrix(features = granges(ilc.integrated),
                                  pwm = pwm_set,
                                  genome = 'mm10',
                                  use.counts = FALSE)

motif.object <- CreateMotifObject(data = motif.matrix,
                                  pwm = pwm_set)

ilc.integrated <- SetAssayData(ilc.integrated,
                               assay = 'ATAC',
                               slot = 'motifs',
                               new.data = motif.object)

ilc.integrated <- RunChromVAR(
  object = ilc.integrated,
  genome = BSgenome.Mmusculus.UCSC.mm10
)


# Save Integration Seurat Object -----------------------------------------------

saveRDS(ilc.integrated,
        file = paste0(out.dir, "pre-transplant_post-transplant_integration_replicate1.rds"))


# Save Session Info ------------------------------------------------------------

writeLines(capture.output(sessionInfo()), "./sessionInfo/integrate_multiome_scRNA_replicate2_sessionInfo.txt")
