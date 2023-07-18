# Integrates 15 human scATAC-seq data sets


library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library(EnsDb.Hsapiens.v86)

plan("multisession", workers = 20)
options(future.globals.maxSize = 50000 * 1024^2) # for 50 Gb RAM


# Reference GRanges for hg38 ---------------------------------------------------
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"


# Output Directory -------------------------------------------------------------
out.dir <- "./Figure_4/"


# Create a common peak set -----------------------------------------------------

# read in each of the peak sets for the human scATAC data
experiments <- c("I","B", "J",
                 "C","D","E",
                 "F","G", "H",
                 "K","L","M",
                 "N","O","P")

# Read in peak Sets
peakSet_files <- list()

for (experiment in experiments) {
  
  df <- read.table(
    file = sprintf("./processed_data/sample%s_peaks.bed", experiment),
    col.names = c("chr", "start", "end") 
  )
  
  peakSet_files <- append(peakSet_files, list(df))
  
}

# convert to genomic ranges
gr_peakSets <- list()

for (peakset in peakSet_files) {
  
  gr <- makeGRangesFromDataFrame(peakset)
  
  gr <- keepStandardChromosomes(gr, pruning.mode = "coarse")
  
  gr_peakSets <- append(gr_peakSets, list(gr))
}

# Create a union set of peaks to quantify in each dataset
combined_peakSet <- reduce(x = c(gr_peakSets[[1]],gr_peakSets[[2]],gr_peakSets[[3]],
                                 gr_peakSets[[4]],gr_peakSets[[5]],gr_peakSets[[6]],
                                 gr_peakSets[[7]],gr_peakSets[[8]],gr_peakSets[[9]],
                                 gr_peakSets[[10]],gr_peakSets[[11]],gr_peakSets[[12]],
                                 gr_peakSets[[13]],gr_peakSets[[14]],gr_peakSets[[15]]
))

# Filter out bad peaks based on length
peakwidths <- width(combined_peakSet)
combined_peakSet <- combined_peakSet[peakwidths  < 10000 & peakwidths > 20]
combined_peakSet


# Create Fragment objects ------------------------------------------------------

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
    cells = rownames(meta_data)
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
                                      genome = 'hg38',
                                      min.cells = 10,
                                      annotation = annotations)
  
  seuratObjectStore <- CreateSeuratObject(chrom_assay, assay = "ATAC", meta.data = meta_data[[i]])
  
  seurat_objects <- append(seurat_objects, list(seuratObjectStore))
  
}

#### Merge Objects ####
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

# Save Merged Seurat Object ----------------------------------------------------
saveRDS(combined, file = paste0(out.dir,"human_scATAC_Merged_SeuratObject.rds"))


# Save Session Info ------------------------------------------------------------
writeLines(capture.output(sessionInfo()), "./sessionInfo/integrate_human_scATAC_sessionInfo.txt")
