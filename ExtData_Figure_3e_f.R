# Generates Extended Data Figure 3e-f


library(ggrepel)
library(tidyverse)
library(Seurat)
library(Signac)


# Output Directory -------------------------------------------------------------
out.dir <- "./Figure_3/"


# Read Replicate 1 Data - Seurat Object ----------------------------------------
ilc.integrated <- readRDS(file = "./Figure_3/inVitro_ILC2_pcILC2_integrated_rep1.rds")


# Identify ILC2 and pcILC2 regulators based on exp condition -------------------

# Identify differential genes
markers_rna <- presto:::wilcoxauc.Seurat(X = ilc.integrated,
                                         group_by = 'orig.ident',
                                         assay = 'data',
                                         seurat_assay = 'SCT')

markers_motifs <- presto:::wilcoxauc.Seurat(X = ilc.integrated,
                                            group_by = 'orig.ident',
                                            assay = 'data',
                                            seurat_assay = 'chromvar')

motif.names <- markers_motifs$feature

colnames(markers_rna) <- paste0("RNA.", colnames(markers_rna))

colnames(markers_motifs) <- paste0("motif.", colnames(markers_motifs))

markers_rna$gene <- markers_rna$RNA.feature

markers_motifs$gene <- ConvertMotifID(ilc.integrated, id = motif.names)


# MA Plot of differential RNA --------------------------------------------------

# Select ex.ILC2 data to plot
data <- markers_rna %>% dplyr::filter(RNA.group == 'pcILC2s') %>% arrange(RNA.logFC)

ma.data <- data %>% mutate(sig.color  = ifelse(-log10(RNA.padj) > -log10(0.05), '#1F9FE1', '#E0601E')) %>%
  arrange(sig.color)

ILC2.pcILC2.RNA.ma.plot <- ggplot(ma.data, aes(x=log2(RNA.avgExpr + 1), y=RNA.logFC,  col = sig.color, label = gene)) + 
  geom_point(alpha = .51) +
  scale_color_identity() +
  geom_text_repel(box.padding = 0.5,
                  size=3,
                  segment.size = 0.25,
                  segment.alpha = 0.5,
                  force_pull=25,
                  max.overlaps=5,
                  min.segment.length= 0.01,
                  data=subset(ma.data, (RNA.feature == "Bcl2" | RNA.feature == "Inpp4b" |
                                          RNA.feature == "Stat1" | RNA.feature == "Il5" |
                                          RNA.feature == "Il13")),
                  colour = 'black') +
  theme_classic() + 
  theme(legend.position="none") +
  geom_hline(yintercept=0, linetype="dashed", color = "black") + 
  ylab("RNA Abundance \n log2(pcILC2 / ILC2)") +
  xlab("RNA Abundance \n log2(Avg. Expression + 1)")

ggsave(ILC2.pcILC2.RNA.ma.plot,
       file = paste0(out.dir, "ExtData_Figure_3e_ILC2_pcILC2_differentialRNA_maPLot.pdf"),
       height = 4.5,
       width = 4.5,
       device = "pdf")


## Save number of differential genes -------------------------------------------
num.diff.genes <- ma.data %>% dplyr::filter(sig.color == "#1F9FE1") %>% nrow()

write("Number of Differential genes:",
      file= paste0(out.dir, "ExtData_Figure_3e_numberDifferentialGenes.txt"))

write(capture.output(num.diff.genes),
      file= paste0(out.dir, "ExtData_Figure_3e_numberDifferentialGenes.txt"),
      append=TRUE)


# Volcano Plot of differential Motif -------------------------------------------

# Select ex.ILC2 data to plot
motif.data <- markers_motifs %>% dplyr::filter(motif.group == 'pcILC2s') %>% arrange(desc(motif.logFC))

motif.ma.data <- motif.data %>% mutate(sig.color  = ifelse(-log10(motif.padj) > -log10(0.05), '#1F9FE1', '#E0601E')) %>%
  arrange(sig.color)


ILC2.pcILC2.motif.volcano.plot <- ggplot(motif.ma.data, aes(x=motif.logFC, y=-log10(motif.padj),  col = sig.color, label = gene)) + 
  geom_point() +
  scale_color_identity() +
  geom_text_repel(box.padding = .5,
                  size=3,
                  segment.size = 0.25,
                  force_pull=1,
                  max.overlaps=500,
                  min.segment.length= 0.01,
                  data=subset(motif.ma.data, (gene == "REL" | gene == "POU5F1" |
                                                gene == "POU2F3" | gene == "IRF7" |
                                                gene == "IRF8" | gene == "IRF2" |
                                                gene == "IRF4")),
                  colour = 'black') +
  theme_classic() + 
  theme(legend.position="none") +
  ylab("-log10(p-value)") +
  xlab("Motif Enrichment \n log2(pcILC2 / ILC2)") +
  xlim(-3,3)

ggsave(ILC2.pcILC2.motif.volcano.plot,
       file = paste0(out.dir, "ExtData_Figure_3f_ILC2_pcILC2_differentialMotif_volcanoPlot.pdf"),
       height = 4.5,
       width = 4.5,
       device = "pdf")

# Save Session Info ------------------------------------------------------------
writeLines(capture.output(sessionInfo()), "./sessionInfo/ExtData_Figure_3e_f_sessionInfo.txt")

