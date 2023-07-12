# Generates Figure 3d


library(tidyverse)
library(ggrepel)


# Output Directory -------------------------------------------------------------
out.dir <- "./Figure_3/"


# Read RT2 Data ----------------------------------------------------------------

rt2.up.data <- read.table("./processed_data/rt2_data_pdf_upreg.txt",
                          header = FALSE,
                          skip = 1)

rt2.down.data <- read.table("./processed_data/rt2_data_pdf_downreg.txt",
                            header = FALSE)

# Set Column names
colnames(rt2.up.data) <- c("Well", "Gene", "FoldChange", "pval", "cat_no")
colnames(rt2.down.data) <- c("Well", "Gene", "FoldChange", "pval", "cat_no")

# log2 transform foldchange upregulated
rt2.up.data$FoldChange <- log2(rt2.up.data$FoldChange)

#Calculate log10 p value  upregulated 
rt2.up.data$neg.log10.pval <- -log10(rt2.up.data$pval)

# log2 transform foldchange downregulated
rt2.down.data$FoldChange <- abs(rt2.down.data$FoldChange)

rt2.down.data$FoldChange <- log2(rt2.down.data$FoldChange)

rt2.down.data$FoldChange <- rt2.down.data$FoldChange * -1

#Calculate log10 p value  down regulated 
rt2.down.data$neg.log10.pval <- -log10(rt2.down.data$pval)

# Concat data frames
rt2.data <- rbind(rt2.up.data,rt2.down.data)

# Add up and down annotation 
rt2.data$sig.status <- with(rt2.data, ifelse(pval < 0.05 & FoldChange > 0,'exILC2',
                                             ifelse(pval <0.05 & FoldChange <0, 'ILC2','nosig')))

rt2.volcano.plot <- ggplot(rt2.data, aes(x = FoldChange, y = neg.log10.pval, colour = sig.status, label = Gene)) +
  geom_point() +
  theme_classic() +
  scale_colour_manual(values = c("#f0972b","#2c5d73",'gray')) +
  theme(legend.position="none") +
  aes(stroke = 0.25) +
  geom_vline(xintercept = 1, linetype="dotted") +
  geom_vline(xintercept = -1, linetype="dotted") +
  ylim(0,7.5) +
  geom_text_repel(box.padding = .5,
                  size=3,
                  segment.size = 0.25,
                  force_pull=1,
                  max.overlaps=50,
                  min.segment.length= 0.01,
                  data=subset(rt2.data, (Gene == "Ifng" | Gene == "Stat1" | Gene == "Fasl" | Gene == "Tbx21" |
                                           Gene == "Il13" | Gene == "Il5" | Gene == "Tnfrsf8")),
                  colour = 'black')

ggsave(rt2.volcano.plot,
       file = paste0(out.dir, "Fig3d_rt2_expression_volcanoPlot.pdf"),
       height = 8,
       width = 5.5 ,
       device = "pdf")


# Save Session Info ------------------------------------------------------------
writeLines(capture.output(sessionInfo()), "./sessionInfo/Figure_3d_sessionInfo.txt")

