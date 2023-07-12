# Identify promoters with differential regions of chromatin accessibility
# ATAC-seq


library(DESeq2)


# Output Directory -------------------------------------------------------------

out.dir <- "./Figure_3/"


# Read in the GBA H3K4me3 count data ------------------------------------------
ATAC.count.matrix <- read.csv(paste0(out.dir,"ILC2_pcILC2_count_matrix.txt"),
                                 sep = '\t',
                                 header = )

# Update header
names(ATAC.count.matrix) = c('geneName', 
                                'ILC1_rep1', 'ILC1_rep2',
                                'ILC2_rep1', 'ILC2_rep2',
                                'ILC2_rep3')

# Set row names as gene names
row.names(ATAC.count.matrix) = ATAC.count.matrix$geneName

# Remove geneName column
ATAC.count.matrix$geneName <- NULL


# Create coldata data frame ----------------------------------------------------
exp_names <- c("ILC1_rep1", "ILC1_rep2", "ILC2_rep1",
               "ILC2_rep2", "ILC2_rep3")

treatment <- c("exp1","exp1","exp2","exp2","exp2")

batch <- c("A", "B", "A", "B", "C")

coldata <- data.frame(treatment, batch,
                      row.names = exp_names,
                      stringsAsFactors = TRUE)


# Identify significantly different H3K4me3 in at least one ILC subtype ---------

dds <- DESeqDataSetFromMatrix(ATAC.count.matrix, coldata, ~treatment)

dds = DESeq(dds, betaPrior=T)

ILC2_pcILC2_results = results(dds, contrast=c("treatment", "exp1", "exp2"))

# Write ILC_LRT_deseq2_results -------------------------------------------------
write.table(ILC2_pcILC2_results,
            file = paste0(out.dir, "pcILC2_ILC2_deseq2_results.txt"),
            row.names=T,
            col.names=NA,
            sep="\t",
            quote = F)


write.table(ILC2_pcILC2_results[ILC2_pcILC2_results$padj < 0.05 & !is.na(ILC2_pcILC2_results$padj),],
            file = paste0(out.dir, "pcILC2_ILC2_deseq2_results_diff_0.05.txt"),
            row.names=T,
            col.names=NA,
            sep="\t",
            quote=F)


# Save Session Info ------------------------------------------------------------
writeLines(capture.output(sessionInfo()), "./sessionInfo/ILC2_pcILC2_DESeq2_sessionInfo.txt")

