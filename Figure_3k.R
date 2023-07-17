# Generates Figure 3k

library(tidyverse)
library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library(EnsDb.Mmusculus.v79)
library(data.table)

plan("multisession", workers = 4)
options(future.globals.maxSize = 50000 * 1024^2) # for 50 Gb RAM


# Reference GRanges for mm10 ---------------------------------------------------
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"


# Output Directory -------------------------------------------------------------
output.dir <- "./Figure_3/"


# Load H5 Raw Data -------------------------------------------------------------

ILC2_data <- Read10X_h5("./processed_data/ilc2_pre-transplant_rep1_filtered_feature_bc_matrix.h5")

pcILC2_data <- Read10X_h5("./processed_data/pcILC2_rep1_filtered_feature_bc_matrix.h5")

post_transplant_data <- Read10X_h5("./processed_data/ILC2_post-transplant_rep1_filtered_feature_bc_matrix.h5")


# Merged Count matrix for RNA-seq data -----------------------------------------

# Extract RNA counts
ILC2_rna_counts = ILC2_data$`Gene Expression`

pcILC2_rna_counts = pcILC2_data$`Gene Expression`

post_transplant_rna_counts = post_transplant_data$`Gene Expression`

# Create Seurat Objects
ILC2 <- CreateSeuratObject(counts = ILC2_rna_counts, project = "ILC2")

pcILC2 <- CreateSeuratObject(counts = pcILC2_rna_counts, project = "pcILC2")

post_transplant <- CreateSeuratObject(counts = post_transplant_rna_counts, project = "post_transplant")


# Create merged RNA count matrix -----------------------------------------------

ilc.merged <- merge(ILC2, y = c(pcILC2, post_transplant), add.cell.ids = c("ILC2", "pcILC2", "post_transplant"), project = "ILCs")

# Mito percentage
ilc.merged[["percent.mt"]] <- PercentageFeatureSet(ilc.merged, pattern = "^mt-")


# Create merged ATAC count matrix ----------------------------------------------

# Create a common peak set 
  ILC2_peaks <- read.table(
    file = "./processed_data/ILC2_atac_peaks_rep1.bed",
    col.names = c("chr", "start", "end") 
  )
  
  pcILC2_peaks <- read.table(
    file = "./processed_data/pcILC2_atac_peaks_rep1.bed",
    col.names = c("chr", "start", "end") 
  )
  
  post_transplant_peaks <- read.table(
    file = "./processed_data/post_transplant_atac_peaks_rep1.bed",
    col.names = c("chr", "start", "end") 
  ) 

  
# All Peak Sets
peakSet_files <- list(ILC2_peaks, pcILC2_peaks, post_transplant_peaks)


# convert to genomic ranges
gr_peakSets <- list()

for (peakset in peakSet_files) {
  
  gr <- makeGRangesFromDataFrame(peakset)
  
  gr <- keepStandardChromosomes(gr, pruning.mode = "coarse")
  
  gr_peakSets <- append(gr_peakSets, list(gr))
}

# Create a union set of peaks to quantify in each dataset
combined_peakSet <- GenomicRanges::reduce(x = c(gr_peakSets[[1]],gr_peakSets[[2]],gr_peakSets[[3]]))

# Filter out bad peaks based on length
peakwidths <- width(combined_peakSet)
combined_peakSet <- combined_peakSet[peakwidths  < 10000 & peakwidths > 20]


# Create Fragment objects ------------------------------------------------------

meta_data <- list()

fragments <- list()

# Barcode metrics
per_barcode_metrics <- c("./processed_data/ILC2_per_barcode_metrics_rep1.csv",
                  "./processed_data/pcILC2_per_barcode_metrics_rep1.csv",
                  "./processed_data/post_transplant_per_barcode_metrics_rep1.csv")
  
# ATAC fragment files 
fragment_files <- c("./processed_data/ilc2_pre-transplant_rep1_atac_fragments.tsv.gz",
                  "./processed_data/pcILC2_rep1_atac_fragments.tsv.gz",
                  "./processed_data/ILC2_post-transplant_rep1_atac_fragments.tsv.gz")

for (i in 1:length(per_barcode_metrics)) {
  
  # load metadata
  metaData <- read.table(
    file = per_barcode_metrics[i],
    stringsAsFactors = FALSE,
    sep = ",",
    header = TRUE,
    row.names = 1
  )[-1, ]
  
  # perform an initial filtering of low count cells
  # Description of column atac_fragments: 
  # https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/pipelines/latest/output/per_barcode_metrics
  metaData <- metaData[metaData$atac_fragments > 750, ]
  
  # create fragment objects
  frags <- CreateFragmentObject(
    path = fragment_files[i],
    cells = rownames(metaData)
  )
  
  # store meta data
  meta_data <- append(meta_data, list(metaData))
  
  # store fragment object
  fragments <- append(fragments, list(frags))
  
}


# Quantify peaks in each data set ----------------------------------------------

countMatricies <- list()

for (i in 1:length(fragments)) {
  
  counts <- FeatureMatrix(
    fragments = fragments[[i]],
    features = combined_peakSet,
    cells = rownames(meta_data[[i]])
  )
  
  countMatricies <- append(countMatricies, list(counts))
  
}


# Create the objects -----------------------------------------------------------

chromatin_assay <- list()

seurat_objects <- list()

for (i in 1:length(countMatricies)) {
  
  chrom_assay <- CreateChromatinAssay(counts = countMatricies[[i]],
                                      fragments = fragments[[i]],
                                      genome = 'mm10',
                                      min.cells = 10,
                                      annotation = annotations)
  
  seuratObjectStore <- CreateSeuratObject(chrom_assay, assay = "ATAC", meta.data = meta_data[[i]])
  
  seurat_objects <- append(seurat_objects, list(seuratObjectStore))
  
  
}


# Merge Objects ----------------------------------------------------------------
ILC2_atac <- seurat_objects[[1]]
pcILC2_atac <- seurat_objects[[2]]
post_transplant_atac <- seurat_objects[[3]]


# add information to identify dataset of origin 
ILC2_atac$dataset <- "ILC2"
pcILC2_atac$dataset <- "pcILC2"
post_transplant_atac$dataset <- "post_transplant"


# merge all datasets, adding a cell ID to make sure cell names are unique
combined_atac <- merge(
  x = ILC2_atac,
  y = list(pcILC2_atac, post_transplant_atac),
  add.cell.ids = c("ILC2", "pcILC2", "post_transplant")
)


# Add ATAC data to merged object -----------------------------------------------

# Shared cells between RNA and ATAC
shared_cells <- intersect(colnames(ilc.merged),colnames(combined_atac))

ilc_merged_match <- ilc.merged[,shared_cells]

combined_atac_match <- combined_atac[,shared_cells]

ilc_merged_match[["ATAC"]] <- combined_atac_match[["ATAC"]]


# QC --------------------------------------------------------------------------- 

ilc_merged_match <- subset(
  x = ilc_merged_match,
  subset = nCount_ATAC < 7e4 &
    nCount_ATAC > 5e3 &
    nCount_RNA < 25000 &
    nCount_RNA > 1000 &
    percent.mt < 20
)


# RNA analysis -----------------------------------------------------------------
DefaultAssay(ilc_merged_match) <- "RNA"
ilc_merged_match <- SCTransform(ilc_merged_match, verbose = FALSE) %>% RunPCA() %>%
  RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')


# ATAC analysis ----------------------------------------------------------------
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(ilc_merged_match) <- "ATAC"

ilc_merged_match <- RunTFIDF(ilc_merged_match)

ilc_merged_match <- FindTopFeatures(ilc_merged_match, min.cutoff = 'q0')

ilc_merged_match <- RunSVD(ilc_merged_match)

ilc_merged_match <- RunUMAP(ilc_merged_match,
                            reduction = 'lsi',
                            dims = 2:50,
                            reduction.name = "umap.atac",
                            reduction.key = "atacUMAP_")


ilc_merged_match <- FindMultiModalNeighbors(ilc_merged_match, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))

ilc_merged_match <- RunUMAP(ilc_merged_match, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")

ilc_merged_match <- FindClusters(ilc_merged_match, graph.name = "wsnn", algorithm = 3, verbose = FALSE)


# Save Merged Seurat Object ----------------------------------------------------
saveRDS(ilc_merged_match, file = paste0(output.dir,"ILC2_pcILC2_post_transplantILC2_norm_seuratObject.rds"))


# Update orig.ident to include ILC1-like cells ---------------------------------

# In vitro integrated
inVitro_ilc_integrated <- readRDS(file = "./Figure_3/inVitro_ILC2_pcILC2_integrated_rep1.rds")

ilc1_like_cells <- inVitro_ilc_integrated@meta.data %>% dplyr::filter(seurat_clusters == 9) %>% rownames()

ILC1.like.cellIDs <- paste0("pcILC2_",str_sub(ilc1_like_cells, end=-3))

ilc_merged_match@meta.data <- within(ilc_merged_match@meta.data, orig.ident[row.names(ilc_merged_match@meta.data) %in% ILC1.like.cellIDs] <- 'ILC1-like')


# Plot UMAP --------------------------------------------------------------------

ilc_merged_match$orig.ident <- factor(ilc_merged_match$orig.ident,
                                      levels = c("ILC2","pcILC2","post_transplant","ILC1-like"))

merged_wnn_plot <- DimPlot(ilc_merged_match,
                           reduction = "wnn.umap",
                           group.by = "orig.ident",
                           label = FALSE,
                           label.size = 2.5,
                           order = c("ILC1-like","pcILC2"),
                           repel = TRUE,
                           cols = c("#72b28b" ,"#2C5D73","#F0972B","purple")) + 
  aes(stroke = 0.35) +
  ggtitle("") +
  theme(legend.position = c(0.6, 0.15)) +
  ylab("UMAP_2") +
  xlab("UMAP_1")
  
ggsave(merged_wnn_plot,
       file = paste0(output.dir,"Fig3k_ILC2_exILC2_post_transplant_integration.pdf"),
       height = 4.5,
       width = 4.5,
       device = "pdf")

# Save Session Info ------------------------------------------------------------
writeLines(capture.output(sessionInfo()), "./sessionInfo/Figure_3k_sessionInfo.txt")
