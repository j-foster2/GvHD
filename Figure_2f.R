# Generates plot for Figure 2f


library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)


# Output Directory -------------------------------------------------------------

out.dir <- "./Figure_2/"


# Read Integrated Seurat Object ------------------------------------------------

ilc.integrated <- readRDS(paste0(out.dir,"pre-transplant_post-transplant_integration_replicate2.rds"))

# Calculate Average ILC1/2/3 LDG signal on integration -------------------------

#Pokrovskii et al. - Fine-subset of nk associated genes
nk.genes <- c('sct_Il12rb2','sct_Klra9','sct_Klra3','sct_Il10rb','sct_Klra1',
              'sct_Ccl4','sct_Ifng','sct_Gzma','sct_Klrc1','sct_Klra7',
              'sct_Ccr2','sct_Il2rb','sct_Tnfsf8','sct_Tgfb1', 'sct_Il21r',
              'sct_Lef1','sct_Egr2','sct_Eomes','sct_Snai1', 'sct_Klrb1c',
              'sct_Ccl3','sct_Ccl4','sct_Klf2','sct_Klri2','sct_Gzma',
              'sct_Ccl5','sct_Il10ra', 'sct_Klf8','sct_Nr4a2','sct_Ifng',
              'sct_Tbx21','sct_Il21r')

nk.genes <- unique(nk.genes)

#Pokrovskii et al. - Fine-subset of ILC1 associated genes
ilc1.genes <- c('sct_Ccl5','sct_Epas1','sct_Mxd4','sct_Il10','sct_Il12rb2',
                'sct_Snai3','sct_Meis3','sct_Ccl4','sct_Mxd1',
                'sct_Tgfb1','sct_Stat3','sct_Foxc2','sct_Gzma','sct_Klrb1c',
                'sct_Ccl3','sct_Ccl4','sct_Klf2','sct_Klri2','sct_Gzma',
                'sct_Ccl5','sct_Il10ra', 'sct_Klf8','sct_Nr4a2','sct_Ifng',
                'sct_Tbx21','sct_Il21r')

ilc1.genes <- unique(ilc1.genes)

#Pokrovskii et al. - Fine-subset of ILC2 associated genes
ilc2.genes <- c('sct_Bmp2','sct_Il13','sct_Il4','sct_Cxcl1','sct_Ifnar1',
                'sct_Il6','sct_Il10ra','sct_Cxcl3','sct_Tnfsf10','sct_Klrg1',
                'sct_Tnfsf14','sct_Il2ra','sct_Bmp7','sct_Areg','sct_E2f1',
                'sct_Ahr','sct_Klf5','sct_Gata3','sct_Il17rb','sct_Il5',
                'sct_Pparg','sct_Il9r','sct_Stat1','sct_Pou2f2')

ilc2.genes <- unique(ilc2.genes)

#Pokrovskii et al. - Fine-subset of ILC3 associated genes
ilc3.genes <- c('sct_Il1r1','sct_Il17a','sct_Cxcr6','sct_Bmp2','sct_Klrb1b',
                'sct_Cd83','sct_Tnfsf11','sct_Il2ra','sct_Il22','sct_Ccr6',
                'sct_Iltifb','sct_Il17f','sct_Il23r','sct_Il17re','sct_Ccr1',
                'sct_Junb','sct_Nfatc2','sct_Sp3','sct_Zfp105','sct_Maff',
                'sct_Rorc','sct_Nr1d1')

ilc3.genes <- unique(ilc3.genes)

innate.lymphoid.genes <- list(nk.genes, ilc1.genes, ilc2.genes, ilc3.genes)

ilc.type <- c("nk","1", "2", "3")

for (j in 1:length(innate.lymphoid.genes)) {
  
  #Extract RNA-seq signal
  gene.signal = FetchData(object = ilc.integrated, vars = innate.lymphoid.genes[[j]])
  #
  #Calculate per cell average of ILC1 genes
  gene.signal.avg = gene.signal %>%
    mutate(ILC.avg = rowMeans(.))
  
  #Add ILC average signal to meta.data data frame
  ilc.integrated@meta.data[ , ncol(ilc.integrated@meta.data) + 1] <- gene.signal.avg$ILC.avg
  
  #Update column name
  colnames(ilc.integrated@meta.data)[ncol(ilc.integrated@meta.data)] <- paste0("ILC", ilc.type[j], ".genes.avg")
  
}

# FIG 2f - Expression of average lineage defining genes ------------------------

# Set order of conditions
ilc.integrated$orig.ident <- factor(ilc.integrated$orig.ident,
                                    levels = c("pre_transplant",'post_transplant'))

# Set order of ILC signature
avg_LDG <- c("ILC2.genes.avg", "ILCnk.genes.avg", "ILC1.genes.avg", "ILC3.genes.avg")

avg.ldg.plot <-  FeaturePlot(ilc.integrated,
                             features = avg_LDG,
                             order = TRUE,
                             cols = c("lightgrey", "#0e3f60"),
                             keep.scale = "all",
                             split.by ="orig.ident",
                             reduction ='umap') + 
  theme(legend.position = "right") & 
  aes(stroke = 0.25)

ggsave(avg.ldg.plot, file = paste0(out.dir, "Fig2f_pre_post_transplant_ILC_Signature_genes.pdf"),
       height = 8, width = 6 , device = "pdf")

# Save Session Info ------------------------------------------------------------

writeLines(capture.output(sessionInfo()), "./sessionInfo/Figure_2f_sessionInfo.txt")
