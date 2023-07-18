# Generates Extended Data Figure 4a 


library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library(EnsDb.Hsapiens.v86)
library(tidyr)
library(ggplot2)
library(dplyr)
library(data.table)


# Output Directory -------------------------------------------------------------
out.dir <- "./Figure_4/"


# Read scATAC data -------------------------------------------------------------

# hILC2
combined.ILC2 <- readRDS(paste0(out.dir,"human_scATAC_Merged_SeuratObject_hILC2Sig.rds"))

# pc-hILC2
combined.pcILC2 <- readRDS(paste0(out.dir,"human_scATAC_Merged_SeuratObject_pc-hILC2Sig.rds"))

# Random
combined.rand <- readRDS(paste0(out.dir,"human_scATAC_Merged_SeuratObject_RandomSig.rds"))


# Order of samples for plotting ------------------------------------------------

sample.order <- c("Healthy1","Healthy2","ILC2",
                  "PreNoGvHD1","PreNoGvHD2","PreNoGvHD3",
                  "PreGvHD1","PreGvHD2","PreGvHD3",
                  "PostNoGvHD1","PostNoGvHD2","PostNoGvHD3",
                  "PostGvHD1","PostGvHD2","PostGvHD3")


# Distribution of ATAC signal at random peaks shared by ILC2 and pcILC2s -------

# Split data by sample 
combined.split <- SplitObject(combined.rand, split.by = "dataset")

avg.ILC2.sig.cell.sample <- list()

for (atac.data in combined.split) {
  
  # Average ATAC signal across ILC2 associated peaks per cell
  avgILC2.sig.cell <- colMeans(atac.data@assays[["ATAC"]]@data)
  
  # Convert matrix to data frame
  avg.ILC2.sig.cell.df <- as.data.frame(avgILC2.sig.cell)
  
  # Annotate ATAC signal with sampleID
  avg.ILC2.sig.cell.df$sample <-  as.character(atac.data$dataset[1])
  
  # Store the Average ATAC signal across ILC2 associated peaks per cell for each sample
  avg.ILC2.sig.cell.sample <- append(avg.ILC2.sig.cell.sample, list(avg.ILC2.sig.cell.df))
  
}

# Average ATAC signal across ILC2 associated peaks per cell for all samples
ILC2.sig <- do.call(rbind, unname(avg.ILC2.sig.cell.sample))

# 5th percentile of all ATAC singal for random - noise floor 
rand.noise.floor <- quantile(ILC2.sig$avgILC2.sig.cell, probs = c(0.10))

# log10 transform that average ATAC signal at ILC2 associated peaks
ILC2.sig$log10sig <- log10(ILC2.sig$avgILC2.sig.cell + 1)

# Set sample factor levels to order violin plot
ILC2.sig$sample <- factor(ILC2.sig$sample, levels = sample.order)

# Violin plot of Avg. ATAC signal at ILC2 associated peaks
rand_open_plot <- ggplot(data=ILC2.sig, aes(x=sample, y=log10sig)) +
  geom_violin() +
  stat_summary(fun=mean, geom="point", size=2, color="black") +
  theme_classic() +
  geom_hline(yintercept=rand.noise.floor, linetype="dashed", color = "orange")

ggsave(rand_open_plot,
       file = paste0(out.dir, "ExtData_Fig4a_human_scATAC_ATAC_Signal_randomOpenChromatin", ".pdf"),
       height = 3.5,
       width = 7.5,
       device = "pdf")


# Select Cells that are above the random noise threshold -----------------------
ILC2.sig.pos <- ILC2.sig %>% dplyr::filter(avgILC2.sig.cell > rand.noise.floor)

ILC2.sig.pos.cells <- rownames(ILC2.sig.pos)


# ILC2 signature distribution --------------------------------------------------
ILC2.split <- SplitObject(combined.ILC2, split.by = "dataset")

avg.ILC2.sig.cell.sample <- list()

for (atac.data in ILC2.split) {
  
  # Average ATAC signal across ILC2 associated peaks per cell
  avgILC2.sig.cell <- colMeans(atac.data@assays[["ATAC"]]@data)
  
  # Convert matrix to data frame
  avg.ILC2.sig.cell.df <- as.data.frame(avgILC2.sig.cell)
  
  # Annotate ATAC signal with sampleID
  avg.ILC2.sig.cell.df$sample <-  as.character(atac.data$dataset[1])
  
  # Store the Average ATAC signal across ILC2 associated peaks per cell for each sample
  avg.ILC2.sig.cell.sample <- append(avg.ILC2.sig.cell.sample, list(avg.ILC2.sig.cell.df))
  
}

# Average ATAC signal across ILC2 associated peaks per cell for all samples
ILC2.sig <- do.call(rbind, unname(avg.ILC2.sig.cell.sample))

# Subset for cells with positive ATAC signal 
ILC2.sig <- ILC2.sig[rownames(ILC2.sig) %in% ILC2.sig.pos.cells, ]

# log10 transform that average ATAC signal at ILC2 associated peaks
ILC2.sig$log10sig <- log10(ILC2.sig$avgILC2.sig.cell + 1)

# Set sample factor levels to order violin plot
ILC2.sig$sample <- factor(ILC2.sig$sample, levels = sample.order)

# Violin plot of Avg. ATAC signal at ILC2 associated peaks
ILC2_sig_plot <- ggplot(data=ILC2.sig, aes(x=sample, y=log10sig)) +
  geom_violin() +
  stat_summary(fun=mean, geom="point", size=2, color="black") +
  theme_classic()

ggsave(ILC2_sig_plot,
       file = paste0(out.dir, "Fig4h_human_scATAC_ILCSignature.pdf"),
       height = 3.5,
       width = 7.5,
       device = "pdf")


# hILC2 - Statistics -----------------------------------------------------------

write("Linear Model results for for Figure 4:\n",
      file= paste0(out.dir, "Fig4_stats.txt"))

ILC2.sig.stats <- ILC2.sig

ILC2.sig.stats$condition1 <- with(ILC2.sig.stats, ifelse(sample == "PreNoGvHD1" | sample == "PreNoGvHD2"| sample == "PreNoGvHD3", 'pre_no',
                                                         ifelse(sample == "PreGvHD1" | sample == "PreGvHD2" |sample == "PreGvHD3", 'pre_yes',
                                                                ifelse(sample == "PostNoGvHD1" | sample == "PostNoGvHD2" |sample == "PostNoGvHD3", 'post_no',
                                                                       ifelse(sample == "PostGvHD1" | sample == "PostGvHD2" |sample == "PostGvHD3", 'post_yes',
                                                                              ifelse(sample == "Healthy1" | sample == "Healthy2", 'healthy','ILC2'))))))

# Labels for comparison across pre- and post-transplant subgroups
ILC2.sig.stats$condition1 <- factor(ILC2.sig.stats$condition1, levels=c("post_no","post_yes","pre_no","pre_yes","healthy","ILC2"))

# Labels for Pre- v. Post-transplant Samples
ILC2.sig.stats$condition2 <- with(ILC2.sig.stats, ifelse(sample == "PreNoGvHD1" | sample == "PreNoGvHD2"| sample == "PreNoGvHD3", 'pre',
                                                         ifelse(sample == "PreGvHD1" | sample == "PreGvHD2" |sample == "PreGvHD3", 'pre',
                                                                ifelse(sample == "PostNoGvHD1" | sample == "PostNoGvHD2" |sample == "PostNoGvHD3", 'post',
                                                                       ifelse(sample == "PostGvHD1" | sample == "PostGvHD2" |sample == "PostGvHD3", 'post',
                                                                              ifelse(sample == "Healthy1" | sample == "Healthy2", 'healthy','ILC2'))))))

ILC2.sig.stats$condition2 <- factor(ILC2.sig.stats$condition2, levels=c("pre","post","healthy","ILC2"))

# Model on Broader Categories (i.e pre- vs. post-transplant samples)

# Save Figure 4 Stats
model <- lm(log10sig ~ condition1, data = ILC2.sig.stats)

write("Figure 4h Regression and statustics:",
      file= paste0(out.dir, "Fig4_stats.txt"),
      append=TRUE)

ILC2_anova <- anova(model)

write("\nFigure 4h ANOVA",
      file= paste0(out.dir, "Fig4_stats.txt"),
      append=TRUE)

write(capture.output(ILC2_anova),
      file= paste0(out.dir, "Fig4_stats.txt"),
      append=TRUE)

ILC2_summary_stats <- summary(model)

write("\nFigure 4h Linear model coefficients and p-vals",
      file= paste0(out.dir, "Fig4_stats.txt"),
      append=TRUE)

write(capture.output(ILC2_summary_stats),
      file= paste0(out.dir, "Fig4_stats.txt"),
      append=TRUE)


# pc-hILC2 signature distribution ----------------------------------------------
pcILC2.split <- SplitObject(combined.pcILC2, split.by = "dataset")

avg.pchILC2.sig.cell.sample <- list()

for (atac.data in pcILC2.split) {
  
  # Average ATAC signal across ILC2 associated peaks per cell
  avg.pchILC2.sig.cell <- colMeans(atac.data@assays[["ATAC"]]@data)
  
  # Convert matrix to data frame
  avg.ILC2.sig.cell.df <- as.data.frame(avg.pchILC2.sig.cell)
  
  # Annotate ATAC signal with sampleID
  avg.ILC2.sig.cell.df$sample <-  as.character(atac.data$dataset[1])
  
  # Store the Average ATAC signal across ILC2 associated peaks per cell for each sample
  avg.pchILC2.sig.cell.sample <- append(avg.pchILC2.sig.cell.sample, list(avg.ILC2.sig.cell.df))
  
}

# Average ATAC signal across ILC2 associated peaks per cell for all samples
pch.ILC2.sig <- do.call(rbind, unname(avg.pchILC2.sig.cell.sample))

# Subset for cells with positive ATAC signal 
pch.ILC2.sig <- pch.ILC2.sig[rownames(pch.ILC2.sig) %in% ILC2.sig.pos.cells, ]

# log10 transform that average ATAC signal at ILC2 associated peaks
pch.ILC2.sig$log10sig <- log10(pch.ILC2.sig$avg.pchILC2.sig.cell + 1)

# Set sample factor levels to order violin plot
pch.ILC2.sig$sample <- factor(pch.ILC2.sig$sample, levels = sample.order)

# Violin plot of Avg. ATAC signal at ILC2 associated peaks
pch.ILC2_sig_plot <- ggplot(data=pch.ILC2.sig, aes(x=sample, y=log10sig)) +
  geom_violin() +
  stat_summary(fun=mean, geom="point", size=2, color="black") +
  theme_classic()

ggsave(pch.ILC2_sig_plot,
       file = paste0(out.dir, "Fig4i_human_scATAC_pc-hILCSignature.pdf"),
       height = 3.5,
       width = 7.5,
       device = "pdf")


# pc-hILC2 - Statistics --------------------------------------------------------
pcILC2.sig.stats <- pch.ILC2.sig

pcILC2.sig.stats$condition1 <- with(pcILC2.sig.stats, ifelse(sample == "PreNoGvHD1" | sample == "PreNoGvHD2"| sample == "PreNoGvHD3", 'pre_no',
                                                         ifelse(sample == "PreGvHD1" | sample == "PreGvHD2" |sample == "PreGvHD3", 'pre_yes',
                                                                ifelse(sample == "PostNoGvHD1" | sample == "PostNoGvHD2" |sample == "PostNoGvHD3", 'post_no',
                                                                       ifelse(sample == "PostGvHD1" | sample == "PostGvHD2" |sample == "PostGvHD3", 'post_yes',
                                                                              ifelse(sample == "Healthy1" | sample == "Healthy2", 'healthy','ILC2'))))))

# Labels for comparison across pre- and post-transplant subgroups
pcILC2.sig.stats$condition1 <- factor(pcILC2.sig.stats$condition1, levels=c("post_no","post_yes","pre_no","pre_yes","healthy","ILC2"))

# Labels for Pre- v. Post-transplant Samples
pcILC2.sig.stats$condition2 <- with(pcILC2.sig.stats, ifelse(sample == "PreNoGvHD1" | sample == "PreNoGvHD2"| sample == "PreNoGvHD3", 'pre',
                                                         ifelse(sample == "PreGvHD1" | sample == "PreGvHD2" |sample == "PreGvHD3", 'pre',
                                                                ifelse(sample == "PostNoGvHD1" | sample == "PostNoGvHD2" |sample == "PostNoGvHD3", 'post',
                                                                       ifelse(sample == "PostGvHD1" | sample == "PostGvHD2" |sample == "PostGvHD3", 'post',
                                                                              ifelse(sample == "Healthy1" | sample == "Healthy2", 'healthy','ILC2'))))))

pcILC2.sig.stats$condition2 <- factor(pcILC2.sig.stats$condition2, levels=c("pre","post","healthy","ILC2"))

# Model on Broader Categories (i.e pre- vs. post-transplant samples)
model <- lm(log10sig ~ condition1, data = pcILC2.sig.stats)


write("Figure 4i Regression and statustics:",
      file= paste0(out.dir, "Fig4_stats.txt"),
      append=TRUE)

pchILC2_anova <- anova(model)

write("\nFigure 4h ANOVA",
      file= paste0(out.dir, "Fig4_stats.txt"),
      append=TRUE)

write(capture.output(pchILC2_anova),
      file= paste0(out.dir, "Fig4_stats.txt"),
      append=TRUE)

pchILC2_summary_stats <- summary(model)

write("\nFigure 4i Linear model coefficients and p-vals",
      file= paste0(out.dir, "Fig4_stats.txt"),
      append=TRUE)

write(capture.output(pchILC2_summary_stats),
      file= paste0(out.dir, "Fig4_stats.txt"),
      append=TRUE)

# Save Session Info ------------------------------------------------------------
writeLines(capture.output(sessionInfo()), "./sessionInfo/Figure4h_i_ExtData_Fig4a_sessionInfo.txt")

