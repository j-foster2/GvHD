# Calculate ATAC signal at peaks in hILC2 and pc-hILC2 union set 


library(Rsubread)


# Output Directory -------------------------------------------------------------

out.dir <- "./Figure_4/"


# Collect ATAC BAM files ----------------------------------------------------

hILC2_bams <- list.files(path = "./processed_data",
                         pattern = glob2rx("hILC*ATAC*bam"),
                         full.names = T)

pc_hILC2_bams <- list.files(path = "./processed_data",
                        pattern = glob2rx("pc-h*bam"),
                        full.names = T)

ATAC_bams <- c(hILC2_bams, pc_hILC2_bams)


# Read in BED of TSSs for all genes --------------------------------------------

fig4b_peak_unionSet <- read.csv("./Figure_4/pc-hILC2_hILC2_union_200bp.bed",
                        sep = '\t',
                        header = T)

# FAIRE-seq/ATAC-seq/ChIP-seq ROE strand doesn't matter
fig4b_peak_unionSet$Strand <- 'na' 

# Column names
names(fig4b_peak_unionSet) <- c('GeneID','Chr', 'Start', 'End', 'Strand')

# Convert to data frame
fig4b_peak_unionSet = as.data.frame(fig4b_peak_unionSet)


# Generate count matrix for peak in hILC2 and pc-hILC2 union set ---------------

ATAC.count.matrix = featureCounts(ATAC_bams,
                                     annot.ext = fig4b_peak_unionSet,
                                     nthreads = 4,
                                  isPairedEnd=TRUE)


# Write count matrix -----------------------------------------------------------

write.table(ATAC.count.matrix$counts,
            file = paste0(out.dir,'hILC2_pc-hILC2_count_matrix.txt'),
            row.names=T,
            col.names=NA,
            sep="\t",
            quote=F)


# Save Session Info ------------------------------------------------------------

writeLines(capture.output(sessionInfo()), "./sessionInfo/generate_hILC2_pc-hILC2_ATAC_count_matrix_sessionInfo.txt")
