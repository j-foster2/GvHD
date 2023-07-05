# Identify promoters with differential TSSs (including 300 bp upstream and 500
# downstream) with differential H3K4me3 signal. Here we performed LRT to identify 
# TSSs with variable signal regardless of change in direction or sample. 
#
# GBA stands for Gury-BenAri et al. 2016 Cell paper (PMID: 27545347)


library(DESeq2)


# Output Directory -------------------------------------------------------------

out.dir <- "./Figure_1/"


# Read in the GBA H3K4me3 count data ------------------------------------------
H3K4me3.count.matrix <- read.csv(paste0(out.dir,"GBA_H3K4me3_count_matrix.txt"),
                             sep = '\t',
                             header = )

# Update header
names(H3K4me3.count.matrix) = c('geneName', 
                                'ILC1_rep1', 'ILC1_rep2',
                                'ILC2_rep1', 'ILC2_rep2',
                                'ILC3_rep1', 'ILC3_rep2')

# Set row names as gene names
row.names(H3K4me3.count.matrix) = H3K4me3.count.matrix$geneName

# Remove geneName column
H3K4me3.count.matrix$geneName <- NULL


# Create coldata data frame ----------------------------------------------------
condition <- c("ILC1", "ILC1", "ILC2", "ILC2", "ILC3", "ILC3")

is_ILC1 <- c("Yes","Yes","No","No","No","No")

is_ILC2 <- c("No", "No", "Yes", "Yes", "No", "No")

is_ILC3 <- c("No", "No", "No", "No", "Yes", "Yes")

exp_names <- c("ILC1_rep1","ILC1_rep2",
               "ILC2_rep1","ILC2_rep2",
               "ILC3_rep1", "ILC3_rep2")

coldata <- data.frame(condition, is_ILC1, is_ILC2, is_ILC3, exp_names,
                      row.names = exp_names,
                      stringsAsFactors = TRUE)


# Identify significantly different H3K4me3 in at least one ILC subtype ---------

dds <- DESeqDataSetFromMatrix(H3K4me3.count.matrix, coldata, ~condition)

dds = DESeq(dds, test="LRT", reduced= ~1)

ILC_results = results(dds)

# Write ILC_LRT_deseq2_results -------------------------------------------------
write.table(ILC_results,
            file = paste0(out.dir, "ILC_LRT_deseq2_results.txt"),
            row.names=T,
            col.names=NA,
            sep="\t",
            quote = F)

# Save Session Info ------------------------------------------------------------

writeLines(capture.output(sessionInfo()), "./sessionInfo/GBA_H3K4me3_DESeq2_sessionInfo.txt")
