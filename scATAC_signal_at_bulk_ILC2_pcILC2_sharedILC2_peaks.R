# Calculate scATAC signal at peaks identified with bulk ATAC data (ILC2 and pcILC2s)


library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library(EnsDb.Hsapiens.v86)
library(tidyr)
library(ggplot2)
library(dplyr)
library(data.table)


plan("multisession", workers = 20)
options(future.globals.maxSize = 75000 * 1024^2) # for 50 Gb RAM


# Output Directory -------------------------------------------------------------
out.dir <- "./Figure_4/"


# Reference GRanges for hg38 ---------------------------------------------------
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"


# Create a common peak set -----------------------------------------------------

# read in each of the peak sets for the human scATAC data
experiments <- c("I","B", "J",
                 "C","D","E",
                 "F","G", "H",
                 "K","L","M",
                 "N","O","P")


# Read ILC2 Signature peaks ----------------------------------------------------

ILC2.sig.bed <- read.table(file = paste0(out.dir,"huILC2_ATAC_sig.bed"),
                           sep = '\t',
                           header = F)

colnames(ILC2.sig.bed) <- c("chr","start","end","strand")

ILC2_combined_peakSet <- makeGRangesFromDataFrame(ILC2.sig.bed)

# Filter out bad peaks based on length
peakwidths <- width(ILC2_combined_peakSet)

ILC2_combined_peakSet <- ILC2_combined_peakSet[peakwidths  < 10000 & peakwidths > 20]


# Read pcILC2 Signature peaks --------------------------------------------------

pcILC2.sig.bed <- read.table(file = paste0(out.dir, "huexILC2_ATAC_sig.bed"),
                             sep = '\t',
                             header = F)

colnames(pcILC2.sig.bed) <- c("chr","start","end","strand")

pcILC2_combined_peakSet <- makeGRangesFromDataFrame(pcILC2.sig.bed)

# Filter out bad peaks based on length
peakwidths <- width(pcILC2_combined_peakSet)

pcILC2_combined_peakSet <- pcILC2_combined_peakSet[peakwidths  < 10000 & peakwidths > 20]


# Read Random Signature peaks --------------------------------------------------

unchanged.sig.bed <- read.table(file = paste0(out.dir, "hILC2_pc-hILC2_unchanged_ATAC_peaks.bed"),
                                sep = '\t',
                                header = F)

# Randomly select peaks from the ILC2 - pcILC2 shared peak set
set.seed(19)

unchanged.sig.bed <- sample_n(unchanged.sig.bed, dim(ILC2.sig.bed)[1])

colnames(unchanged.sig.bed) <- c("chr","start","end")

unchanged_combined_peakSet <- makeGRangesFromDataFrame(unchanged.sig.bed)

# Filter out bad peaks based on length
peakwidths <- width(unchanged_combined_peakSet)

unchanged_combined_peakSet <- unchanged_combined_peakSet[peakwidths  < 10000 & peakwidths > 20]


# Calculate ATAC signal at peaks within each set -------------------------------

# ATAC signature objects
peak_sets <- list(ILC2_combined_peakSet,
               pcILC2_combined_peakSet,
               unchanged_combined_peakSet)

# RDS output files 
out_files <- list("human_scATAC_Merged_SeuratObject_hILC2Sig.rds",
               "human_scATAC_Merged_SeuratObject_pc-hILC2Sig.rds",
               "human_scATAC_Merged_SeuratObject_RandomSig.rds")

for (j in 1:length(out_files)) {
  
  # Create Fragment objects ----------------------------------------------------

  meta_data <- list()

  fragments <- list()

  for (experiment in experiments) {

    # load metadata
    metaData <- read.table(
      file = sprintf("./processed_data/sample%s_singlecell.csv", experiment),
      stringsAsFactors = FALSE,
      sep = ",",
      header = TRUE,
      row.names = 1
    )[-1, ]

    # perform an initial filtering of low count cells
    metaData <- metaData[metaData$passed_filters > 750, ]

    # create fragment objects
    frags <- CreateFragmentObject(
      path = sprintf("./processed_data/sample%s_fragments.tsv.gz", experiment),
      cells = rownames(metaData)
    )

    # store meta data
    meta_data <- append(meta_data, list(metaData))

    # store fragment object
    fragments <- append(fragments, list(frags))

  }


  # Quantify peaks in each dataset ---------------------------------------------

  countMatricies <- list()

  for (i in 1:length(fragments)) {

    counts <- FeatureMatrix(
      fragments = fragments[[i]],
      features = peak_sets[[j]],
      cells = rownames(meta_data[[i]])
    )

    countMatricies <- append(countMatricies, list(counts))

  }


  # Create the objects ---------------------------------------------------------
  chromatin_assay <- list()

  seurat_objects <- list()

  for (i in 1:length(countMatricies)) {

    chrom_assay <- CreateChromatinAssay(counts = countMatricies[[i]],
                                        fragments = fragments[[i]],
                                        genome = 'hg38',
                                        min.cells = 10,
                                        annotation = annotations)

    seuratObjectStore <- CreateSeuratObject(chrom_assay, assay = "ATAC", meta.data = meta_data[[i]])

    seurat_objects <- append(seurat_objects, list(seuratObjectStore))

  }


  # Merge Objects --------------------------------------------------------------
  ILC2 <- seurat_objects[[1]]
  Healthy_Donor1 <- seurat_objects[[2]]
  Healthy_Donor2 <- seurat_objects[[3]]
  PreTreatment_noGVHD1 <- seurat_objects[[4]]
  PreTreatment_noGVHD2 <- seurat_objects[[5]]
  PreTreatment_noGVHD3 <- seurat_objects[[6]]
  PreTreatment_GVHD1 <- seurat_objects[[7]]
  PreTreatment_GVHD2 <- seurat_objects[[8]]
  PreTreatment_GVHD3 <- seurat_objects[[9]]
  PostTreatment_noGVHD1 <- seurat_objects[[10]]
  PostTreatment_noGVHD2 <- seurat_objects[[11]]
  PostTreatment_noGVHD3 <- seurat_objects[[12]]
  PostTreatment_GVHD1 <- seurat_objects[[13]]
  PostTreatment_GVHD2 <- seurat_objects[[14]]
  PostTreatment_GVHD3 <- seurat_objects[[15]]

  # add information to identify dataset of origin
  ILC2$dataset <- "ILC2"
  Healthy_Donor1$dataset <- "Healthy1"
  Healthy_Donor2$dataset <- "Healthy2"
  PreTreatment_noGVHD1$dataset <- "PreNoGvHD1"
  PreTreatment_noGVHD2$dataset <- "PreNoGvHD2"
  PreTreatment_noGVHD3$dataset <- "PreNoGvHD3"
  PreTreatment_GVHD1$dataset <- "PreGvHD1"
  PreTreatment_GVHD2$dataset <- "PreGvHD2"
  PreTreatment_GVHD3$dataset <- "PreGvHD3"
  PostTreatment_noGVHD1$dataset <- "PostNoGvHD1"
  PostTreatment_noGVHD2$dataset <- "PostNoGvHD2"
  PostTreatment_noGVHD3$dataset <- "PostNoGvHD3"
  PostTreatment_GVHD1$dataset <- "PostGvHD1"
  PostTreatment_GVHD2$dataset <- "PostGvHD2"
  PostTreatment_GVHD3$dataset <- "PostGvHD3"

  # merge all datasets, adding a cell ID to make sure cell names are unique
  combined <- merge(
    x = ILC2,
    y = list(Healthy_Donor1, Healthy_Donor2,
             PreTreatment_noGVHD1, PreTreatment_noGVHD2, PreTreatment_noGVHD3,
             PreTreatment_GVHD1, PreTreatment_GVHD2, PreTreatment_GVHD3,
             PostTreatment_noGVHD1, PostTreatment_noGVHD2, PostTreatment_noGVHD3,
             PostTreatment_GVHD1, PostTreatment_GVHD2, PostTreatment_GVHD3),
    add.cell.ids = c("ILC2", "Healthy1", "Healthy2",
                     "PreNo1","PreNo2","PreNo3",
                     "PreYes1","PreYes2","PreYes3",
                     "PostNo1", "PostNo2","PostNo3",
                     "PostYes1","PostYes2","PostYes3")
  )
  combined[["ATAC"]]


  # Normalization, dim reduction, UMAP -----------------------------------------
  combined <- RunTFIDF(combined)

  combined <- FindTopFeatures(combined, min.cutoff = 20)

  combined <- RunSVD(combined)

  combined <- RunUMAP(combined, dims = 2:50, reduction = 'lsi')

  
  # Save Merged Seurat Object --------------------------------------------------
  saveRDS(combined, file = paste0(out.dir,out_files[j]))
  
}


# Save Session Info ------------------------------------------------------------
writeLines(capture.output(sessionInfo()), "./sessionInfo/scATAC_signal_at_bulk_ILC2_pcILC2_randomILC2_peaks_sessionInfo.txt")

