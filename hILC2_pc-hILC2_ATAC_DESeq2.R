# Identify promoters with differential regions of chromatin accessibility
# ATAC-seq


library(DESeq2)


# Output Directory -------------------------------------------------------------

out.dir <- "./Figure_4/"


# Read in the GBA H3K4me3 count data ------------------------------------------
ATAC.count.matrix <- read.csv(paste0(out.dir,"hILC2_pc-hILC2_count_matrix.txt"),
                              sep = '\t',
                              header = )

# Update header
names(ATAC.count.matrix) = c('geneName', 
                             'hILC2_rep1', 'hILC2_rep2',
                             'pc-hILC2_rep1', 'pc-hILC2_rep2')

# Set row names as gene names
row.names(ATAC.count.matrix) = ATAC.count.matrix$geneName

# Remove geneName column
ATAC.count.matrix$geneName <- NULL


# Create coldata data frame ----------------------------------------------------
exp_names <- c("hILC2_rep1", "hILC2_rep2",
               "pc-hILC2_rep1", "pc-hILC2_rep2")

treatment <- c("exp1","exp1","exp2","exp2")

batch <- c("A", "B", "A", "B")

coldata <- data.frame(treatment, batch,
                      row.names = exp_names,
                      stringsAsFactors = TRUE)


# Identify significantly different ATAC ---------

dds <- DESeqDataSetFromMatrix(ATAC.count.matrix, coldata, ~treatment)

dds = DESeq(dds, betaPrior=T)

hILC2_pc_hILC2_results = results(dds, contrast=c("treatment", "exp1", "exp2"))

# Write deseq2_results -------------------------------------------------
write.table(hILC2_pc_hILC2_results,
            file = paste0(out.dir, "hILC2_pc-hILC2_deseq2_results.txt"),
            row.names=T,
            col.names=NA,
            sep="\t",
            quote = F)


write.table(hILC2_pc_hILC2_results[hILC2_pc_hILC2_results$padj < 0.05 & !is.na(hILC2_pc_hILC2_results$padj),],
            file = paste0(out.dir, "hILC2_pc-hILC2_deseq2_results_diff_0.05.txt"),
            row.names=T,
            col.names=NA,
            sep="\t",
            quote=F)


# Save Session Info ------------------------------------------------------------
writeLines(capture.output(sessionInfo()), "./sessionInfo/hILC2_pc-hILC2_ATAC_DESeq2_sessionInfo.txt")

