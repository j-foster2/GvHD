# Generates plot for Figure 2d
#
# Calculate overlap of gene sets identified through Figure 2d differential gene
# expression analysis.


library(UpSetR)


# Output Directory -------------------------------------------------------------

out.dir <- "./Figure_2/"


# Read list of gene sets from Fig. 2d ------------------------------------------

diff_gene_list <- readRDS(paste0(out.dir,"Fig2d_gene_sets.rds"))

# Extended Data Fig 2d ---------------------------------------------------------

# Copy the gene sets associated with each cluster
gene_sets_exp <- diff_gene_list

# Give cluster name to each gene set
names(gene_sets_exp) <- c("Pre1","Pre2","Pre3", "Post1", "Post2","Post3")

# Upset Plot
pdf(paste0(out.dir, "ExtData_Fig2d_geneSets_overlap.pdf"),
    onefile=FALSE)

upset(fromList(gene_sets_exp), nsets = 6,
      sets = c("Post3", "Post2","Post1","Pre3","Pre2","Pre1"),
      order.by = c("freq"),
      keep.order = TRUE)

dev.off()

# Save Session Info ------------------------------------------------------------

writeLines(capture.output(sessionInfo()), "./sessionInfo/ExtData_Figure_2d_sessionInfo.txt")
