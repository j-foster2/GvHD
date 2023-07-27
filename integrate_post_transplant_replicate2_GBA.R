# Integrates post-transplant single-nucleus multiome data with previously
# published MARS-seq data of CD127+ cell population from the small intestinal
# lamina propria of healthy mice (PMID: 27545347). This integration is based on
# the RNA data only.
#
# Memory Requirement: 100g


library(Seurat)
library(patchwork)
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(Signac)
library(EnsDb.Mmusculus.v79)
library(GenomicAlignments)
library(chromVAR)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(BSgenome.Mmusculus.UCSC.mm10)
library(BiocParallel)

# Output Directory -------------------------------------------------------------
out.dir <- "./Figure_2/"


# Load mm10 genome data --------------------------------------------------------
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)

seqlevelsStyle(annotations) <- 'UCSC'

genome(annotations) <- "mm10"


# Load motif database ----------------------------------------------------------
pwm_set <- getMatrixSet(x = JASPAR2020,
                        opts = list(collection = "CORE",
                                    tax_group = 'vertebrates',
                                    all_versions = FALSE))


# Load Multiome data - HDF5 format ---------------------------------------------

# Load Post-transplant
post_transplant_data <- Read10X_h5("./processed_data/ILC2_post-transplant_rep2_filtered_feature_bc_matrix.h5")


# Gury-BenAri data -------------------------------------------------------------

# Load in Gury-BenAri UMI Table
Gury.BenAri.UMI = read.table(file = paste0("./GBA_Processed_Data/GSE85152_umitab_noQuote.txt"),
                             sep = '\t',
                             row.names = 1,
                             quote = "",
                             header = T)

# Isolate single geneID for each row
single.gene.annotaiton = str_split_fixed(rownames(Gury.BenAri.UMI), ";", 2)[,1]

rownames(Gury.BenAri.UMI) = single.gene.annotaiton

# Create Seurat object
ILCs = CreateSeuratObject(counts = Gury.BenAri.UMI, project = "Gury_BenAri")

ILCs <- subset(
  x = ILCs,
  subset = 
    nCount_RNA > 200)

# Post-transplant data ---------------------------------------------------------

# Extract RNA and ATAC count data
post_transplant_rna_counts = post_transplant_data$`Gene Expression`

post_transplant_atac_counts <- post_transplant_data$Peaks

# Create Seurat object for post_transplant RNA count data
post_transplant = CreateSeuratObject(counts = post_transplant_rna_counts, project = "post_transplant")

# Add  Post-transplant ATAC-seq data (standard chromosomes)
grange.counts <- StringToGRanges(rownames(post_transplant_atac_counts), sep = c(":", "-"))

grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)

post_transplant_atac_counts <- post_transplant_atac_counts[as.vector(grange.use), ]

# ATAC fragment file
frag.file <- "./processed_data/ILC2_post-transplant_rep2_atac_fragments.tsv.gz"

# Create ATAC assay object
chrom_assay <- CreateChromatinAssay(
  counts = post_transplant_atac_counts,
  sep = c(":", "-"),
  genome = 'mm10',
  fragments = frag.file,
  min.cells = 10,
  annotation = annotations
)

# Add post-transplant ATAC assay object to seurat object
post_transplant[["ATAC"]] <- chrom_assay


# Post transplant  QC ----------------------------------------------------------

# Calculate mitochondrial reads for each cell
post_transplant[["percent.mt"]] <- PercentageFeatureSet(post_transplant, pattern = "^mt-")

# QC Filter
post_transplant <- subset(
  x = post_transplant,
  subset = nCount_ATAC < 7e4 &
    nCount_ATAC > 5e3 &
    nCount_RNA < 25000 &
    nCount_RNA > 1000 &
    percent.mt < 20)


# ATAC analysis ----------------------------------------------------------------

DefaultAssay(post_transplant) <- "ATAC"

# Term-frequency inverse-document-frequency Normalization
# Explanation: https://stuartlab.org/signac/articles/pbmc_vignette.html
post_transplant <- RunTFIDF(post_transplant)

post_transplant <- FindTopFeatures(post_transplant, min.cutoff = 'q0')

post_transplant <- RunSVD(post_transplant)

post_transplant <- RunUMAP(post_transplant,
                           reduction = 'lsi',
                           dims = 2:50,
                           reduction.name = "umap.atac",
                           reduction.key = "atacUMAP_")

# motif ID
DefaultAssay(post_transplant) <- "ATAC"

motif.matrix <- CreateMotifMatrix(features = granges(post_transplant),
                                  pwm = pwm_set,
                                  genome = 'mm10',
                                  use.counts = FALSE)

motif.object <- CreateMotifObject(data = motif.matrix,
                                  pwm = pwm_set)

post_transplant <- SetAssayData(post_transplant,
                                assay = 'ATAC',
                                slot = 'motifs',
                                new.data = motif.object)

register(MulticoreParam(20))
post_transplant <- RunChromVAR(
  object = post_transplant,
  genome = BSgenome.Mmusculus.UCSC.mm10
)


# Save Post-transplant Seurat Object ------------------------------------------------
saveRDS(post_transplant, file = paste0(out.dir,"Post-transplant_SeuratObject_replicate2.rds"))


# Integrate Post-transplant and GBA samples by RNA (Seurat Workflow) -----------

DefaultAssay(post_transplant) <- "RNA"

DefaultAssay(ILCs) <- "RNA"

#List of Seurat Objects
ILC.list = list(ILCs,post_transplant)


ILC.list = lapply(X = ILC.list, FUN = function(x) {
  x <- SCTransform(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features = SelectIntegrationFeatures(object.list = ILC.list)

immune.anchors = FindIntegrationAnchors(object.list = ILC.list,
                                        anchor.features = features)

ilc.combined = IntegrateData(anchorset = immune.anchors)


DefaultAssay(ilc.combined) <- "integrated"

# Seurat Integration workflow
ilc.combined = ScaleData(ilc.combined, verbose = FALSE)

ilc.combined = RunPCA(ilc.combined, npcs = 30, verbose = FALSE)

ilc.combined = RunUMAP(ilc.combined, reduction = "pca", dims = 1:30)

ilc.combined = FindNeighbors(ilc.combined, reduction = "pca", dims = 1:30)

ilc.combined = FindClusters(ilc.combined, resolution = 0.5)


# Save integrated Seurat Object ------------------------------------------------
saveRDS(ilc.combined, file = paste0(out.dir,"Post-transplant_GBA_integration.rds"))


# Save Session Info ------------------------------------------------------------
writeLines(capture.output(sessionInfo()), "./sessionInfo/integrate_post_transplant_replicate2_GBA_sessionInfo.txt")

