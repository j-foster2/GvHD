# Generate Figure 1b Plot

library(tidyr)
library(dplyr)
library(ggplot2)


# Output Directory -------------------------------------------------------------

out.dir <- "./Figure_1/"

# Read data --------------------------------------------------------------------
data.path = "./Figure_1/"

rna_k4me3_data <- read.table(file = paste0(data.path, "bruce_rna_laurie_k4me3_data.txt"),
                                sep = '\t',
                                header = T)
#drop Na
rna_k4me3_data <- rna_k4me3_data %>% drop_na()


# Calculate correlation between of RNA and K4me3 -------------------------------
spearman_corr <- cor(rna_k4me3_data$ILC2_RNA_avg, rna_k4me3_data$ILC2_K4me3_avg, method = "spearman")

corr.label <- paste0("r = ",round(spearman_corr, digits = 2))

hex.code <- "#3286CC"


# Plot RNA v. K4me3 ------------------------------------------------------------
rna_k4me3_data_plot <- ggplot(log2(rna_k4me3_data), aes(x=ILC2_RNA_avg, y=ILC2_K4me3_avg)) + 
  geom_point(shape=20, color=hex.code, size = 1, alpha = 0.30 ) +
  theme_classic() +
  annotate("text", x = -8.5, y = 5, label = corr.label,, size = 1)

# Save Plot
ggsave(rna_k4me3_data_plot, file = paste0(out.dir, "Bruce_RNA_Laurie_H3K4me3_correlation", ".pdf"),
       height = 6, width = 6 , device = "pdf")

# Save Session Info ------------------------------------------------------------

writeLines(capture.output(sessionInfo()), "./sessionInfo/Figure_1b_sessionInfo.txt")

