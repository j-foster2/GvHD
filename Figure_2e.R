# Generates plot for Figure 2e


library(ggplot2)


# Output Directory -------------------------------------------------------------

out.dir <- "./Figure_2/"


# Read Enrichr Data associated with post-transplant gene sets ------------------

enrichr_data <- read.table(file = paste0(out.dir,"post_transplant_enrichr.txt"),
                           sep = '\t',
                           quote = "",
                           header = T)

enrichr_data$cluster <- as.factor(enrichr_data$cluster)

enrichr_data$cluster <- factor(enrichr_data$cluster, levels=c('Post3',
                                                              'Post2',
                                                              'Post1'))


# Dot Plot ---------------------------------------------------------------------
enrichr_plot <- ggplot(enrichr_data, aes(x=reorder(Term,-value,FUN=mean),
                                         y=cluster,
                                         size=value,
                                         color=overlap )) + 
  geom_point() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) + 
  ylab("") + 
  xlab("") +
  guides(size = guide_legend(title = "-log10(padj)"))

ggsave(enrichr_plot, file = paste0(out.dir,"Fig2e_post_transplant_geneSet_Enrichr_MouseGeneAtlas.pdf"),
       height = 7, width = 4 , device = "pdf")


# Save Session Info ------------------------------------------------------------

writeLines(capture.output(sessionInfo()), "./sessionInfo/Figure_2e_sessionInfo.txt")
