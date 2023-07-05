# Calculate H3K4me3 signal at TSS (including 300 bp upstream and 500 downstream)
#
# GBA stands for Gury-BenAri et al. 2016 Cell paper (PMID: 27545347)

library(Rsubread)


# Output Directory -------------------------------------------------------------

out.dir <- "./Figure_1/"


# Collect H3K4me3 BAM files ----------------------------------------------------

K4me3_bams <- list.files(path = "./processed_data",
                         pattern = glob2rx("blfILC**GBA**bam"),
                         full.names = T)


# Read in BED of TSSs for all genes --------------------------------------------

refSeq_TSSs <- read.csv(paste0(out.dir,"refseq_promotor.bed"),
                        sep = '\t',
                        header = F)

# FAIRE-seq/ATAC-seq/ChIP-seq ROE strand doesn't matter
refSeq_TSSs$Strand <- 'na' 

# Column names
names(refSeq_TSSs) <- c('Chr', 'Start', 'End', 'GeneID', 'Strand')

# Reorder Columns
refSeq_TSSs <- refSeq_TSSs[, c('GeneID', 'Chr', 'Start', 'End', 'Strand')]

# Convert to data frame
refSeq_TSSs = as.data.frame(refSeq_TSSs)


# Generate count matrix of H3K4me3 reads at gene TSSs --------------------------

H3K4me3.count.matrix = featureCounts(K4me3_bams,
                                       annot.ext = refSeq_TSSs,
                                       nthreads = 4)


# Write GBA_K4me3 counts at TSSs -----------------------------------------------

write.table(H3K4me3.count.matrix$counts,
            file = paste0(out.dir,'GBA_H3K4me3_count_matrix.txt'),
            row.names=T,
            col.names=NA,
            sep="\t",
            quote=F)


# Save Session Info ------------------------------------------------------------

writeLines(capture.output(sessionInfo()), "./sessionInfo/generate_GBA_H3K4me3_count_matrix_sessionInfo.txt")
