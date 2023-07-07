# Integrates all multiome data (2 replicates of pre-transplant and 2 replicates
# of post-transplant cells)


library(Seurat)
library(patchwork)
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggrepel)
library(gprofiler2)
library(forcats)

library(Signac)
library(EnsDb.Mmusculus.v79)
library(GenomicAlignments)


#Motif Identification libraries and data
library(chromVAR)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(BSgenome.Mmusculus.UCSC.mm10)

pwm_set <- getMatrixSet(x = JASPAR2020,
                        opts = list(collection = "CORE",
                                    tax_group = 'vertebrates',
                                    all_versions = FALSE))

# Output Directory -------------------------------------------------------------
out.dir <- "./Figure_2/"


# Reference GRanges for mm10 ---------------------------------------------------
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"


# Raw Data ---------------------------------------------------------------------

pre_transplant_rep1 <- Read10X_h5("./processed_data/ilc2_pre-transplant_rep1_filtered_feature_bc_matrix.h5")

pre_transplant_rep2 <- Read10X_h5("./processed_data/ilc2_pre-transplant_rep2_filtered_feature_bc_matrix.h5")

post_transplant_rep1 <- Read10X_h5("./processed_data/ILC2_post-transplant_rep1_filtered_feature_bc_matrix.h5")

post_transplant_rep2 <- Read10X_h5("./processed_data/ILC2_post-transplant_rep2_filtered_feature_bc_matrix.h5")

pre_datasets  <-  list(pre_transplant_rep1, pre_transplant_rep2)

post_datasets  <-  list(post_transplant_rep1, post_transplant_rep2)


# Select Biological Replicate --------------------------------------------------

exp.rep = c("rep1", "rep2")

i <- 1

# Read ILC2 data Replicate 1 - (including peak identification) -----------------

#Extract RNA and ATAC data
pre.rna_counts_rep1 = pre_datasets[[i]]$`Gene Expression`

pre.atac_counts_rep1 <- pre_datasets[[i]]$Peaks

#Create Seurat object for ILC2 data
pre_rep1 = CreateSeuratObject(counts = pre.rna_counts_rep1, project = "pre_rep1")

# Now add in the ATAC-seq data - we'll only use peaks in standard chromosomes
grange.counts <- StringToGRanges(rownames(pre.atac_counts_rep1), sep = c(":", "-"))

grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)

pre.atac_counts_rep1 <- pre.atac_counts_rep1[as.vector(grange.use), ]

# Set fragment file
if (i == 1){
  frag.file <- "./processed_data/ilc2_pre-transplant_rep1_atac_fragments.tsv.gz"
  
} else{
  
  frag.file <- "./processed_data/ilc2_pre-transplant_rep2_atac_fragments.tsv.gz"  
  
}

chrom_assay <- CreateChromatinAssay(
  counts = pre.atac_counts_rep1,
  sep = c(":", "-"),
  genome = 'mm10',
  fragments = frag.file,
  min.cells = 10,
  annotation = annotations
)

pre_rep1[["ATAC"]] <- chrom_assay


# Read GFP data Replicate 1 - (including peak identification) ------------------

#Extract RNA and ATAC data
post.rna_counts_rep1 = post_datasets[[i]]$`Gene Expression`

post.atac_counts_rep1 <- post_datasets[[i]]$Peaks

#Create Seurat object for GFP data
post_rep1 = CreateSeuratObject(counts = post.rna_counts_rep1, project = "post_rep1")

# Now add in the ATAC-seq data - we'll only use peaks in standard chromosomes
grange.counts <- StringToGRanges(rownames(post.atac_counts_rep1), sep = c(":", "-"))

grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)

post.atac_counts_rep1 <- post.atac_counts_rep1[as.vector(grange.use), ]

if (i == 1){
  gfp.frag.file <- "./processed_data/ILC2_post-transplant_rep1_atac_fragments.tsv.gz"
  
} else{
  
  gfp.frag.file <- "./processed_data/ILC2_post-transplant_rep2_atac_fragments.tsv.gz"  
  
}

chrom_assay <- CreateChromatinAssay(
  counts = post.atac_counts_rep1,
  sep = c(":", "-"),
  genome = 'mm10',
  fragments = gfp.frag.file,
  min.cells = 10,
  annotation = annotations
)

post_rep1[["ATAC"]] <- chrom_assay


# Set Biological Replicate -----------------------------------------------------
i <- 2


# Read Pre-transplant data Replicate 2 - (including peak identification) -------

#Extract RNA and ATAC data
pre.rna_counts_rep2 = pre_datasets[[i]]$`Gene Expression`

pre.atac_counts_rep2 <- pre_datasets[[i]]$Peaks

#Create Seurat object for ILC2 data
pre_rep2 = CreateSeuratObject(counts = pre.rna_counts_rep2, project = "pre_rep2")

# Now add in the ATAC-seq data - we'll only use peaks in standard chromosomes
grange.counts <- StringToGRanges(rownames(pre.atac_counts_rep2), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
pre.atac_counts_rep2 <- pre.atac_counts_rep2[as.vector(grange.use), ]

# Set fragment file
if (i == 1){
  frag.file <- "./processed_data/ilc2_pre-transplant_rep1_atac_fragments.tsv.gz"
  
} else{
  
  frag.file <- "./processed_data/ilc2_pre-transplant_rep2_atac_fragments.tsv.gz"  
  
}

chrom_assay <- CreateChromatinAssay(
  counts = pre.atac_counts_rep2,
  sep = c(":", "-"),
  genome = 'mm10',
  fragments = frag.file,
  min.cells = 10,
  annotation = annotations
)

pre_rep2[["ATAC"]] <- chrom_assay


# Read GFP data Replicate 2 - (including peak identification) ------------------

#Extract RNA and ATAC data
post.rna_counts_rep2 = post_datasets[[i]]$`Gene Expression`

post.atac_counts_rep2 <- post_datasets[[i]]$Peaks

#Create Seurat object for GFP data
post_rep2 = CreateSeuratObject(counts = post.rna_counts_rep2, project = "post_rep2")

# Now add in the ATAC-seq data - we'll only use peaks in standard chromosomes
grange.counts <- StringToGRanges(rownames(post.atac_counts_rep2), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
post.atac_counts_rep2 <- post.atac_counts_rep2[as.vector(grange.use), ]

if (i == 1){
  gfp.frag.file <- "./processed_data/ILC2_post-transplant_rep1_atac_fragments.tsv.gz"
  
} else{
  
  gfp.frag.file <- "./processed_data/ILC2_post-transplant_rep2_atac_fragments.tsv.gz"  
  
}

chrom_assay <- CreateChromatinAssay(
  counts = post.atac_counts_rep2,
  sep = c(":", "-"),
  genome = 'mm10',
  fragments = gfp.frag.file,
  min.cells = 10,
  annotation = annotations
)

post_rep2[["ATAC"]] <- chrom_assay


# QC Filtering -----------------------------------------------------------------

pre_rep1[["percent.mt"]] <- PercentageFeatureSet(pre_rep1, pattern = "^mt-")
pre_rep2[["percent.mt"]] <- PercentageFeatureSet(pre_rep2, pattern = "^mt-")
post_rep1[["percent.mt"]] <- PercentageFeatureSet(post_rep1, pattern = "^mt-")
post_rep2[["percent.mt"]] <- PercentageFeatureSet(post_rep2, pattern = "^mt-")

pre_rep1 <- subset(
  x = pre_rep1,
  subset = nCount_ATAC < 7e4 &
    nCount_ATAC > 5e3 &
    nCount_RNA < 25000 &
    nCount_RNA > 1000 &
    percent.mt < 20
)

pre_rep2 <- subset(
  x = pre_rep2,
  subset = nCount_ATAC < 7e4 &
    nCount_ATAC > 5e3 &
    nCount_RNA < 25000 &
    nCount_RNA > 1000 &
    percent.mt < 20
)

post_rep1 <- subset(
  x = post_rep1,
  subset = nCount_ATAC < 7e4 &
    nCount_ATAC > 5e3 &
    nCount_RNA < 25000 &
    nCount_RNA > 1000 &
    percent.mt < 20
)

post_rep2 <- subset(
  x = post_rep2,
  subset = nCount_ATAC < 7e4 &
    nCount_ATAC > 5e3 &
    nCount_RNA < 25000 &
    nCount_RNA > 1000 &
    percent.mt < 20
)

# Integration using RNA Data ------------------------------------------------

DefaultAssay(pre_rep1) <- "RNA"
DefaultAssay(pre_rep2) <- "RNA"

DefaultAssay(post_rep1) <- "RNA"
DefaultAssay(post_rep2) <- "RNA"

# List of Seurat Objects
ILC.list = list(pre_rep1,pre_rep2,post_rep1,post_rep2)

# normalize and identify variable features for each data set independently
ILC.list = lapply(X = ILC.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features = SelectIntegrationFeatures(object.list = ILC.list)

immune.anchors = FindIntegrationAnchors(object.list = ILC.list,
                                        anchor.features = features)

ilc.integrated = IntegrateData(anchorset = immune.anchors)

# Scale -> PCA -> UMAP Integrated data
DefaultAssay(ilc.integrated) <- "integrated"

# Scale integrated Data
ilc.integrated = ScaleData(ilc.integrated, verbose = FALSE)

# PCA on integrated data
ilc.integrated = RunPCA(ilc.integrated,
                        npcs = 30,
                        verbose = FALSE,
                        reduction.name = "pcaIntegrated")

# UMAP on integrated data
ilc.integrated = RunUMAP(ilc.integrated,
                         reduction = "pcaIntegrated",
                         dims = 1:30)

# 
ilc.integrated = FindNeighbors(ilc.integrated,
                               reduction = "pcaIntegrated",
                               dims = 1:30)

# Clustering
ilc.integrated = FindClusters(ilc.integrated,
                              resolution = 0.5)


# Save ILC RNA Integrated Seurat Object ----------------------------------------
saveRDS(ilc.integrated, file = paste0(out.dir,"Pre_Post_transplant_multiome_integration.rds"))


# Save Session Info ------------------------------------------------------------
writeLines(capture.output(sessionInfo()), "./sessionInfo/integrate_all_multiome_data_sessionInfo.txt")
