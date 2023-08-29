# Integration of in vitro mouse ILC2 and pcILC2 data


library(dplyr)
library(Seurat)
library(Signac)
library(patchwork)
library(tidyr)
library(stringr)
library(forcats)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
library(GenomicAlignments)
library(chromVAR)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(BiocParallel)


# Output Directory -------------------------------------------------------------
out.dir <- "./Figure_3/"


# Load Motif Database ----------------------------------------------------------
pwm_set <- getMatrixSet(x = JASPAR2020,
                        opts = list(collection = "CORE",
                                    tax_group = 'vertebrates',
                                    all_versions = FALSE))


# mm10 Genome information ------------------------------------------------------
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"


# Load Raw Data ----------------------------------------------------------------
ILC2.data <- Read10X_h5("./processed_data/ilc2_pre-transplant_rep1_filtered_feature_bc_matrix.h5")

pcILC2.data <- Read10X_h5("./processed_data/pcILC2_rep1_filtered_feature_bc_matrix.h5")


# Read ILC2 data - (including peak identification) -----------------------------

#Extract RNA and ATAC data
ILC2.rna_counts = ILC2.data$`Gene Expression`

ILC2.atac_counts <- ILC2.data$Peaks

# Create Seurat object for ILC2 data
ILC2 = CreateSeuratObject(counts = ILC2.rna_counts, project = "ILC2s")

# Add ATAC-seq data - standard chromosomes
grange.counts <- StringToGRanges(rownames(ILC2.atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
ILC2.atac_counts <- ILC2.atac_counts[as.vector(grange.use), ]

# ATAC Fragment File
frag.file <- "./processed_data/ilc2_pre-transplant_rep1_atac_fragments.tsv.gz"

chrom_assay <- CreateChromatinAssay(
  counts = ILC2.atac_counts,
  sep = c(":", "-"),
  genome = 'mm10',
  fragments = frag.file,
  min.cells = 10,
  annotation = annotations
)

ILC2[["ATAC"]] <- chrom_assay


# Read pcILC2 data - (including peak identification) ---------------------------

# Extract RNA and ATAC data
pcILC2.rna_counts = pcILC2.data$`Gene Expression`

pcILC2.atac_counts <- pcILC2.data$Peaks

# Create Seurat object for ILC2 data
pcILC2 = CreateSeuratObject(counts = pcILC2.rna_counts, project = "pcILC2s")

# Add ATAC-seq data - standard chromosomes
grange.counts <- StringToGRanges(rownames(pcILC2.atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
pcILC2.atac_counts <- pcILC2.atac_counts[as.vector(grange.use), ]

# ATAC fragment file
pcILC2.frag.file <- "./processed_data/pcILC2_rep1_atac_fragments.tsv.gz"

chrom_assay <- CreateChromatinAssay(
  counts = pcILC2.atac_counts,
  sep = c(":", "-"),
  genome = 'mm10',
  fragments = pcILC2.frag.file,
  min.cells = 10,
  annotation = annotations
)

pcILC2[["ATAC"]] <- chrom_assay


# QC Filtering -----------------------------------------------------------------

ILC2[["percent.mt"]] <- PercentageFeatureSet(ILC2, pattern = "^mt-")

pcILC2[["percent.mt"]] <- PercentageFeatureSet(pcILC2, pattern = "^mt-")

ILC2 <- subset(
  x = ILC2,
  subset = nCount_ATAC < 7e4 &
    nCount_ATAC > 5e3 &
    nCount_RNA < 25000 &
    nCount_RNA > 1000 &
    percent.mt < 20
)

pcILC2 <- subset(
  x = pcILC2,
  subset = nCount_ATAC < 7e4 &
    nCount_ATAC > 5e3 &
    nCount_RNA < 25000 &
    nCount_RNA > 1000 &
    percent.mt < 20
)


# Integration (Seurat Workflow) ------------------------------------------------

DefaultAssay(ILC2) <- "RNA"

DefaultAssay(pcILC2) <- "RNA"

ILC.list = list(ILC2, pcILC2)

ILC.list = lapply(X = ILC.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features = SelectIntegrationFeatures(object.list = ILC.list)

immune.anchors = FindIntegrationAnchors(object.list = ILC.list,
                                        anchor.features = features)

ilc.integrated = IntegrateData(anchorset = immune.anchors)

DefaultAssay(ilc.integrated) <- "integrated"

ilc.integrated = ScaleData(ilc.integrated, verbose = FALSE)

ilc.integrated = RunPCA(ilc.integrated, npcs = 30, verbose = FALSE)

ilc.integrated = RunUMAP(ilc.integrated, reduction = "pca", dims = 1:30)

ilc.integrated = FindNeighbors(ilc.integrated, reduction = "pca", dims = 1:30)

ilc.integrated = FindClusters(ilc.integrated, resolution = 0.8)

DimPlot(ilc.integrated)


# RNA analysis -----------------------------------------------------------------

DefaultAssay(ilc.integrated) <- "RNA"

ilc.integrated <- SCTransform(ilc.integrated, verbose = FALSE) %>%
  RunPCA() %>%
  RunUMAP(dims = 1:50,
          reduction.name = 'umap.rna',
          reduction.key = 'rnaUMAP_')


# ATAC analysis ----------------------------------------------------------------

DefaultAssay(ilc.integrated) <- "ATAC"

ilc.integrated <- RunTFIDF(ilc.integrated)

ilc.integrated <- FindTopFeatures(ilc.integrated, min.cutoff = 'q0')

ilc.integrated <- RunSVD(ilc.integrated)

ilc.integrated <- RunUMAP(ilc.integrated,
                          reduction = 'lsi',
                          dims = 2:50,
                          reduction.name = "umap.atac",
                          reduction.key = "atacUMAP_")

# motif identification
DefaultAssay(ilc.integrated) <- "ATAC"

motif.matrix <- CreateMotifMatrix(features = granges(ilc.integrated),
                                  pwm = pwm_set,
                                  genome = 'mm10',
                                  use.counts = FALSE)

motif.object <- CreateMotifObject(data = motif.matrix, pwm = pwm_set)

ilc.integrated <- SetAssayData(ilc.integrated,
                               assay = 'ATAC',
                               slot = 'motifs',
                               new.data = motif.object)

register(MulticoreParam(20))

ilc.integrated <- RunChromVAR(
  object = ilc.integrated,
  genome = BSgenome.Mmusculus.UCSC.mm10
)

# Save Integrated Seurat Object ----------------------------------------
saveRDS(ilc.integrated, file = paste0(out.dir,"inVitro_ILC2_pcILC2_integrated_rep1.rds"))


# Save Session Info ------------------------------------------------------------
writeLines(capture.output(sessionInfo()), "./sessionInfo/inVitro_ILC2_pcILC2_integrated_rep1_R_sessionInfo.txt")

