# Calculate ATAC signal at union set of ATAC peaks identified in mouse ILC2
# and pcILC2 cells

library(Rsubread)


# Output Directory -------------------------------------------------------------

out.dir <- "./Figure_3/"


# Collect H3K4me3 BAM files ----------------------------------------------------

ATAC_bams <- list.files(path = "./processed_data",
                         pattern = glob2rx("mILC**ATAC**bam"),
                         full.names = T)


# Read in BED of TSSs for all genes --------------------------------------------

fig3e_peak_unionSet <- read.csv("./Figure_3/pcILC2_ILC2_union_200bp.bed",
                        sep = '\t',
                        header = T)

# FAIRE-seq/ATAC-seq/ChIP-seq ROE strand doesn't matter
fig3e_peak_unionSet$Strand <- 'na' 

# Column names
names(fig3e_peak_unionSet) <- c('GeneID','Chr', 'Start', 'End', 'Strand')

# Convert to data frame
fig3e_peak_unionSet = as.data.frame(fig3e_peak_unionSet)


# Generate count matrix of H3K4me3 reads at gene TSSs --------------------------

ATAC.count.matrix = featureCounts(ATAC_bams,
                                     annot.ext = fig3e_peak_unionSet,
                                     nthreads = 4,
                                  isPairedEnd=TRUE)


# Write GBA_K4me3 counts at TSSs -----------------------------------------------

write.table(ATAC.count.matrix$counts,
            file = paste0(out.dir,'ILC2_pcILC2_count_matrix.txt'),
            row.names=T,
            col.names=NA,
            sep="\t",
            quote=F)


# Save Session Info ------------------------------------------------------------

writeLines(capture.output(sessionInfo()), "./sessionInfo/generate_ILC2_pcILC2_count_matrix_sessionInfo.txt")
