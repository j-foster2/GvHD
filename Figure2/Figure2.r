#### Figure 2 Code ####
# Joey Foster
# Davis Lab
#######################


library(JASPAR2020)
library(TFBSTools)
library(Seurat)
library(EnsDb.Mmusculus.v79)
library(GenomicAlignments) #Works
library(Signac)
library(patchwork)
library(tidyr)
library(dplyr)
library(stringr) #Works
library(ggplot2)
library(ggrepel)
library(gprofiler2)
library(forcats)
library(clustree)
library(gplots)
library(UpSetR)
library(chromVAR)
library(motifmatchr)
library(BSgenome.Mmusculus.UCSC.mm10) #Works
library(circlize)
library(ComplexHeatmap)


#### Load Motif Database ####
pwm_set <- getMatrixSet(x = JASPAR2020,
                        opts = list(collection = "CORE",
                                    tax_group = 'vertebrates',
                                    all_versions = FALSE))


#### Reference GRanges for mm10 ####
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "mm10"


#### Replicate Choice ####
i <- 2 # Data shown in Fig 2b-d, 2f-m
# i<- 2 # Extended Data Fig 2g-l


#### Output Directory ####
out.dir <- "/proj/dllab/jfoster/serody_project/GitHub/Figure2/"


#### Load Multiome data - HDF5 format ####

#Load Pre-transplant and post-transplant data 
ILC2.data.rep1 <- Read10X_h5("/proj/dllab/jfoster/serody_project/results/tenX_scMultiome/ILC2_rep1/ILC2_rep1/outs/filtered_feature_bc_matrix.h5")

ILC2.data.rep2 <- Read10X_h5("/proj/dllab/jfoster/serody_project/results/tenX_scMultiome/ILC2_rep2/ILC2_rep2/outs/filtered_feature_bc_matrix.h5")

gfp.data.rep1 <- Read10X_h5("/proj/dllab/jfoster/serody_project/results/ex-ILC2_multiome/2021-12-02/exILC2_GVHD_sample1/outs/filtered_feature_bc_matrix.h5")

gfp.data.rep2 <- Read10X_h5("/proj/dllab/jfoster/serody_project/results/ex-ILC2_multiome_rep2/2021-12-09/exILC2_GVHD_sample2/outs/filtered_feature_bc_matrix.h5")


# Gernate lists of samples - used to toggle between replicates 
ILC2.datasets  <-  list(ILC2.data.rep1, ILC2.data.rep2)

gfp.datasets  <-  list(gfp.data.rep1, gfp.data.rep2)


# Track sample being analyzed
exp.rep = c("rep1", "rep2")


#### Read Pre-transplant ILC2 data - (including peak identification) ####

# Extract RNA and ATAC data
ILC2.rna_counts = ILC2.datasets[[i]]$`Gene Expression`

ILC2.atac_counts <- ILC2.datasets[[i]]$Peaks

# Create Seurat object for pre-transplant data
ILC2 = CreateSeuratObject(counts = ILC2.rna_counts, project = "ILC2s")

# Add ATAC data to seurat object
grange.counts <- StringToGRanges(rownames(ILC2.atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
ILC2.atac_counts <- ILC2.atac_counts[as.vector(grange.use), ]

if (i == 1){
  frag.file <- "/proj/dllab/jfoster/serody_project/results/tenX_scMultiome/ILC2_rep1/ILC2_rep1/outs/atac_fragments.tsv.gz"
  
} else{
  
  frag.file <- "/proj/dllab/jfoster/serody_project/results/tenX_scMultiome/ILC2_rep2/ILC2_rep2/outs/atac_fragments.tsv.gz"  
  
}

chrom_assay <- CreateChromatinAssay(
  counts = ILC2.atac_counts,
  sep = c(":", "-"),
  genome = 'mm10',
  fragments = frag.file,
  min.cells = 10,
  annotation = annotations
)

ILC2[["ATAC"]] <- chrom_assay

#### Read Post-transplant ILC2 data - (including peak identification) ####

# Extract RNA and ATAC data
gfp.rna_counts = gfp.datasets[[i]]$`Gene Expression`

gfp.atac_counts <- gfp.datasets[[i]]$Peaks

# Create Seurat object for GFP data
gfp = CreateSeuratObject(counts = gfp.rna_counts, project = "GFP")

# Add ATAC data to Seurat Object
grange.counts <- StringToGRanges(rownames(gfp.atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
gfp.atac_counts <- gfp.atac_counts[as.vector(grange.use), ]

if (i == 1){
  gfp.frag.file <- "/proj/dllab/jfoster/serody_project/results/ex-ILC2_multiome/2021-12-02/exILC2_GVHD_sample1/outs/atac_fragments.tsv.gz"
  
} else{
  
  gfp.frag.file <- "/proj/dllab/jfoster/serody_project/results/ex-ILC2_multiome_rep2/2021-12-09/exILC2_GVHD_sample2/outs/atac_fragments.tsv.gz"  
  
}

#TODO - update chrom_assay object with gfp name
chrom_assay <- CreateChromatinAssay(
  counts = gfp.atac_counts,
  sep = c(":", "-"),
  genome = 'mm10',
  fragments = gfp.frag.file,
  min.cells = 10,
  annotation = annotations
)

gfp[["ATAC"]] <- chrom_assay


#### QC Filtering ####

ILC2[["percent.mt"]] <- PercentageFeatureSet(ILC2, pattern = "^mt-")
# ex.ILC2[["percent.mt"]] <- PercentageFeatureSet(ex.ILC2, pattern = "^mt-")
gfp[["percent.mt"]] <- PercentageFeatureSet(gfp, pattern = "^mt-")

# Filter pre-transplant cells
ILC2 <- subset(
  x = ILC2,
  subset = nCount_ATAC < 7e4 &
    nCount_ATAC > 5e3 &
    nCount_RNA < 25000 &
    nCount_RNA > 1000 &
    percent.mt < 20
)

# Filter post-transplant cells
gfp <- subset(
  x = gfp,
  subset = nCount_ATAC < 7e4 &
    nCount_ATAC > 5e3 &
    nCount_RNA < 25000 &
    nCount_RNA > 1000 &
    percent.mt < 20
)

#### Integrate pre- and post-transplant data ####

DefaultAssay(ILC2) <- "RNA"
DefaultAssay(gfp) <- "RNA"

# List of Seurat Objects
ILC.list = list(ILC2, gfp)

# Normalization and identification of variable features
ILC.list = lapply(X = ILC.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# Identify features that are variable across datasets 
features = SelectIntegrationFeatures(object.list = ILC.list)

immune.anchors = FindIntegrationAnchors(object.list = ILC.list,
                                        anchor.features = features)

# Integrate data using RNA data
ilc.integrated = IntegrateData(anchorset = immune.anchors)

# Run Seurat workflow for Integrated data
DefaultAssay(ilc.integrated) <- "integrated"

# TODO scale data - describe mroe 
ilc.integrated = ScaleData(ilc.integrated, verbose = FALSE)

# Dimension reduction with PCA
ilc.integrated = RunPCA(ilc.integrated, npcs = 30, verbose = FALSE)

# Dimension reduction and visualization of data with UMAP
ilc.integrated = RunUMAP(ilc.integrated,
                         reduction = "pca",
                         dims = 1:30)

# Compute nearest neighbors
ilc.integrated = FindNeighbors(ilc.integrated, reduction = "pca", dims = 1:30)

# Cluster identification at several resolutions
ilc.integrated = FindClusters(ilc.integrated, resolution = c(0, 0.05, 0.1, 0.2, 
                                                             0.3, 0.4, 0.5, 
                                                             0.6, 0.7, 0.8, 
                                                             0.9, 1.0))


#### Evaluate impact of altering resolution for clustering ####

ilc.clusttree <- clustree(ilc.integrated,
                          prefix = "integrated_snn_res.")

ilc.clusttree

# Store the clusttree that shows there are a handful of cells that
# are always in a separate cluster. This holds true even at low 
# resolution values.
ilc.integrated@misc <- append(ilc.integrated@misc, list(ilc.clusttree))

names(ilc.integrated@misc) <- c("clusterOptimizationPlot")

#View plot
ilc.integrated@misc$clusterOptimizationPlot

#### Drop two cells that cluster on their own at multiple resolution values ####

cells.to.remove <- c(row.names(ilc.integrated@meta.data %>% dplyr::filter(ilc.integrated@meta.data['integrated_snn_res.0.05'] == 2)))

ilc.integrated <- ilc.integrated[,!colnames(ilc.integrated) %in% cells.to.remove]


#### RNA analysis (Seurat workflow) ####
DefaultAssay(ilc.integrated) <- "RNA"

ilc.integrated <- SCTransform(ilc.integrated, verbose = FALSE) %>%
  RunPCA() %>%
  RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')


#### ATAC analysis (Seurat workflow) ####

DefaultAssay(ilc.integrated) <- "ATAC"

# Term-frequency inverse-document-frequency Normalization
# Explanation: https://stuartlab.org/signac/articles/pbmc_vignette.html
ilc.integrated <- RunTFIDF(ilc.integrated)

ilc.integrated <- FindTopFeatures(ilc.integrated, min.cutoff = 'q0')

ilc.integrated <- RunSVD(ilc.integrated)

ilc.integrated <- RunUMAP(ilc.integrated, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

# Identification of motifs
DefaultAssay(ilc.integrated) <- "ATAC"

motif.matrix <- CreateMotifMatrix(features = granges(ilc.integrated), pwm = pwm_set, genome = 'mm10', use.counts = FALSE)
motif.object <- CreateMotifObject(data = motif.matrix, pwm = pwm_set)
ilc.integrated <- SetAssayData(ilc.integrated, assay = 'ATAC', slot = 'motifs', new.data = motif.object)

ilc.integrated <- RunChromVAR(
  object = ilc.integrated,
  genome = BSgenome.Mmusculus.UCSC.mm10
)

#### Save ILC RNA Integrated Seurat Object ####

saveRDS(ilc.integrated,
        file = paste0(out.dir, "ILC2_GFP_rna_integrated_",exp.rep[i],"_v4.3", ".rds"))


#### Read ILC RNA Integrated Seurat Object ####

ilc.integrated <- readRDS(paste0(out.dir,
                                 "ILC2_GFP_rna_integrated_",exp.rep[i],"_v4.3", ".rds"))


#### Calculate Average ILC1/2/3 LDG signal on integration ####

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


#### Data visualizations ####


###### Counts of cells in different groups #####

#Number of post-transplant cells that are in pre-transplant clusters
num_post_in_pre <- ilc.integrated@meta.data %>% 
  dplyr::filter(orig.ident == "GFP" &
                  (integrated_snn_res.0.4 == 0 | integrated_snn_res.0.4 == 1 |integrated_snn_res.0.4 == 2)) %>% 
  nrow()

num_post <- ilc.integrated@meta.data %>% 
  dplyr::filter(orig.ident == "GFP" ) %>% nrow()

perct_post_in_pre <- num_post_in_pre / num_post

perct_post_in_pre

#Number of pre-transplant cells that are in post-transplant clusters
num_pre_in_post <- ilc.integrated@meta.data %>% 
  dplyr::filter(orig.ident == "ILC2s" &
                  (integrated_snn_res.0.4 == 3 | integrated_snn_res.0.4 == 4 |integrated_snn_res.0.4 == 5)) %>% 
  nrow()

num_pre <- ilc.integrated@meta.data %>% 
  dplyr::filter(orig.ident == "ILC2s" ) %>% nrow()

perct_pre_in_post <- num_pre_in_post / num_pre

perct_pre_in_post


###### FIG. 2b - UMAP annotated with sample ID ######

ilc.integrated.umap <- DimPlot(ilc.integrated,
                               reduction = "umap",
                               group.by = "orig.ident",
                               label = FALSE,
                               label.size = 2.5,
                               repel = TRUE,
                               cols = c("#2C5D73", "#7cd37b")) + 
  aes(stroke = 0.25) +
  ggtitle("") + 
  theme(legend.position = c(0.01, 0.15))

ilc.integrated.umap

ggsave(ilc.integrated.umap, file = paste0(out.dir,"pre-post_transplant_RNA_integration_UMAP_sampleAnno_noLabels_",exp.rep[i], ".pdf"),
       height = 4.5, width = 4.5 , device = "pdf")


####### FIG. 2c - UMAP -clust. res.0.4 - pre v. post clustering - Split by sample #######
# Website for colors = https://medialab.github.io/iwanthue/
ilc.integrated$orig.ident <- factor(ilc.integrated$orig.ident, levels = c("ILC2s","ex.ILC2s",'GFP'))

ilc.integrated.umap.split <- DimPlot(ilc.integrated,
                                     reduction = "umap",
                                     group.by = 'integrated_snn_res.0.4',
                                     split.by = "orig.ident",
                                     label = FALSE,
                                     label.size = 2.5,
                                     repel = TRUE,
                                     cols = c("#005aa2",
                                              "#62c6c5",
                                              "#a77f48",
                                              "#238654",
                                              "#88B570",
                                              "#56813e")) + aes(stroke = 0.25) +
  ggtitle("") 

ilc.integrated.umap.split

ggsave(ilc.integrated.umap.split, file = paste0(out.dir,"pre-post_transplant_RNA_integration_res-0.4_UMAP-Split_",exp.rep[i], ".pdf"),
       height = 4, width = 8, device = "pdf")


##### FIG 2f - Expression of average lineage defining genes #####

ilc.integrated$orig.ident <- factor(ilc.integrated$orig.ident, levels = c("ILC2s",'GFP'))

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

avg.ldg.plot

ggsave(avg.ldg.plot, file = paste0(out.dir, "pre_post_transplant.ILCSigs_",exp.rep[i], ".pdf"),
       height = 8, width = 6 , device = "pdf")

##### FIG 2g - Distribution of Average LDG across cells - Violin Plot ####

avg.ldg.violin <- VlnPlot(object = ilc.integrated,
        features = avg_LDG,
        cols = c("#2C5D73", "#7cd37b"), 
        split.by ="orig.ident",
        group.by = "orig.ident",
        same.y.lims = TRUE,
        ncol = 1,
        pt.size = 0)

avg.ldg.violin

ggsave(avg.ldg.violin, file = paste0(out.dir, "pre_post_transplant.ILCSigs_violin_",exp.rep[i], ".pdf"),
       height = 10, width =2.5 , device = "pdf")

# Mann Whitney U test

#Pre v. Post - ILC2 Signature 
pre_ilc2 <- ilc.integrated@meta.data %>% dplyr::filter(orig.ident == "ILC2s") %>% dplyr::select("ILC2.genes.avg")

post_ilc2 <- ilc.integrated@meta.data %>% dplyr::filter(orig.ident == "GFP") %>% dplyr::select("ILC2.genes.avg")

wilcox.test(pre_ilc2$ILC2.genes.avg, post_ilc2$ILC2.genes.avg)

#Pre v. Post - NK Signature 
pre_nk <- ilc.integrated@meta.data %>% dplyr::filter(orig.ident == "ILC2s") %>% dplyr::select("ILCnk.genes.avg")

post_nk <- ilc.integrated@meta.data %>% dplyr::filter(orig.ident == "GFP") %>% dplyr::select("ILCnk.genes.avg")

wilcox.test(pre_nk$ILCnk.genes.avg, post_nk$ILCnk.genes.avg)

#Pre v. Post - ILC1 Signature 
pre_ilc1 <- ilc.integrated@meta.data %>% dplyr::filter(orig.ident == "ILC2s") %>% dplyr::select("ILC1.genes.avg")

post_ilc1 <- ilc.integrated@meta.data %>% dplyr::filter(orig.ident == "GFP") %>% dplyr::select("ILC1.genes.avg")

wilcox.test(pre_ilc1$ILC1.genes.avg, post_ilc1$ILC1.genes.avg)

#Pre v. Post - ILC3 Signature 
pre_ilc3 <- ilc.integrated@meta.data %>% dplyr::filter(orig.ident == "ILC2s") %>% dplyr::select("ILC3.genes.avg")

post_ilc3 <- ilc.integrated@meta.data %>% dplyr::filter(orig.ident == "GFP") %>% dplyr::select("ILC3.genes.avg")

wilcox.test(pre_ilc3$ILC3.genes.avg ,
            post_ilc3$ILC3.genes.avg,
            alternative = "two.sided")

##### FIG. 2d - Heatmap - Average Gene Expression #####

###### Differential gene expression analysis of clusters ######

Idents(object = ilc.integrated) <- "integrated_snn_res.0.4"

# Store differential genes for each cluster
diff_gene_list <- list()

# pre-transplant differential genes relative to post-transplant cells
pre.transplant.clusters <- c(0,2,1)

for (pre.clust in pre.transplant.clusters) {
  
  # Identify differential genes
  diff.genes <- FindMarkers(ilc.integrated,
                            ident.1 = pre.clust,
                            ident.2 = c(3,4,5),
                            min.pct = 0.25,
                            assay = "SCT")  
  
  # Select significant genes with FC > 0
  top.clusterX.genes <- diff.genes %>%
    dplyr::filter(p_val_adj < 0.05 & avg_log2FC > 0) %>%
    arrange(desc(avg_log2FC))
  
  diff_gene_list <- append(diff_gene_list,list(row.names(top.clusterX.genes)))
  
  # Write text file of differential genes for each cluster
  write.table(top.clusterX.genes, 
              paste0(out.dir, "ILC2_transplant_differential_genes_cluster_",pre.clust,".txt"),
              sep = '\t',
              row.names = T,
              quote = F)
  
}

# post-transplant differential genes relative to all pre-transplant cells
post.transplant.clusters <- c(3,4,5)

for (post.clust in post.transplant.clusters) {
  
  # Identify differential genes
  diff.genes <- FindMarkers(ilc.integrated,
                            ident.1 = post.clust,
                            ident.2 = c(0,1,2),
                            min.pct = 0.25,
                            assay = "SCT")  
  
  # Select significant genes with FC > 0
  top.clusterX.genes <- diff.genes %>%
    dplyr::filter(p_val_adj < 0.05 & avg_log2FC > 0) %>%
    arrange(desc(avg_log2FC))
  
  diff_gene_list <- append(diff_gene_list, list(row.names(top.clusterX.genes)))
  
  # Write txt file of differential genes for each cluster
  write.table(top.clusterX.genes,paste0(out.dir, "ILC2_transplant_differential_genes_cluster_",post.clust,".txt"),
              sep = '\t',
              row.names = T,
              quote = F)
}

###### Generate Union Set of differential genes  ######

# Clusters to cycle through
cluster.names <- levels(ilc.integrated)

cluster.names

#Need vector of cell IDs for each cluster 
cluster.cellIDs <- list()

#List of matricies of average matricies 
avg_gene_exp <- list()

# List of avg gene dataframes
gene_set_avg <- list()

for (gene_set in diff_gene_list) {
  
  for (cluster in cluster.names) {
    
    # Isolate the cellIDs for each cluster
    cell_IDs <- row.names(ilc.integrated@meta.data %>% dplyr::filter(integrated_snn_res.0.4 == cluster))
    
    # Normalized RNA data
    rna_data <- ilc.integrated@assays[["SCT"]]@data
    
    # Submatrix of normalized RNA signal
    sub_rna_data <- rna_data[gene_set, cell_IDs]
    
    # Vector of avg gene expression
    sub_rna_avg <- rowMeans(sub_rna_data)
    
    # Store avg gene values (across cells) for each cluster
    avg_gene_exp <- append(avg_gene_exp, list(sub_rna_avg))
    
  }
  
  #Concat avgerage gene vectors into dataframe for each gene set
  gene_set_avg_df <-as.data.frame(do.call(cbind, avg_gene_exp))
  
  # Save the 
  gene_set_avg <- append(gene_set_avg, list(gene_set_avg_df))
  
  # Reset the list holding the average score for each gene set
  avg_gene_exp <- list()
  
}

# Highlight pre-transplant shared genes in across pre groups
pre_shared_genes <- intersect(intersect(rownames(gene_set_avg[[1]]),rownames(gene_set_avg[[2]])),rownames(gene_set_avg[[3]]))

# Highlight pre-transplant shared genes in across pre groups
post_shared_genes <- intersect(intersect(rownames(gene_set_avg[[4]]),rownames(gene_set_avg[[5]])),rownames(gene_set_avg[[6]]))

# Copy gene names to column (needed to identified duplicates)
gene_set_avg_keepNames <- list()
for (df in gene_set_avg) {
  
  df$geneName <- rownames(df)
  
  gene_set_avg_keepNames <- append(gene_set_avg_keepNames, list(df))
  
}

# Full Matrix of genes 
diff_gene_clust_data <- do.call(rbind, gene_set_avg_keepNames)

# pre-transplant shared genes boolean of shared or not
pre_shared_genes <- c("pre_not_shared","pre_shared")[1 + (diff_gene_clust_data[,'geneName'] %in% pre_shared_genes)]

# post-transplant shared genes boolean of shared or not
post_shared_genes <- c("post_not_shared","post_shared")[1 + (diff_gene_clust_data[,'geneName'] %in% post_shared_genes)]

# Drop geneName column 
diff_gene_clust_data <- diff_gene_clust_data[-c(7)]

###### Complex Heatmap ######
# Mean Center the rows 
df <- t(scale(t(data.matrix(diff_gene_clust_data)), scale = FALSE))

#Save avg. gene expression dataFrame
# saveRDS(df, file = paste0(out.dir, "geneExpression_Heatmap_data_frame_",exp.rep[i], ".rds"))

# Read avg. gene expression dataFrame 
df <- readRDS(paste0(out.dir, "geneExpression_Heatmap_data_frame_rep2.rds"))

# Remove mean attribute from matrix - Not necessary 
attr(df, "scaled:center") <- NULL

cluster_names = c("Pre1","Pre2","Pre3",
                  "Post1","Post2","Post3")

row_split = rep(cluster_names, times = lengths(diff_gene_list))

#order row split with factor
row_split <- factor(row_split, levels=cluster_names)

pdf(paste0(out.dir, "pre_post_cluster_exp_",exp.rep[i], ".pdf"),
    width=4,
    height=4.25)
  
Heatmap(df,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        row_names_gp = gpar(fontsize = 0),
        name = "Norm.\nGene Exp.",
        left_annotation = HeatmapAnnotation(clusters = row_split,
                                            col = list(clusters = c("Pre1" = "#005ba3",
                                                               "Pre2" = "#a77f48",
                                                               "Pre3" = "#63c7c5",
                                                               "Post1" = "#238755",
                                                               "Post2" = "#88b571",
                                                               "Post3" = "#57813d")),
                                            which = 'row'),
        right_annotation = HeatmapAnnotation(Pre_Intersect = pre_shared_genes,
                                             Post_Intersect = post_shared_genes,
                                             col = list(Pre_Intersect = c("pre_shared" = "#000000",
                                                                     "pre_not_shared" = "#FFFFFF"),
                                                        Post_Intersect = c("post_shared" = "#000000",
                                                                       "post_not_shared" = "#FFFFFF")),
                                            which = 'row',
                                            show_legend = FALSE),
        
        row_split = row_split,
        row_title = NULL)


dev.off()

# Single Cell resolution Heatmap 

#list to vector of differential genes
diff_gene_vec <- unlist(diff_gene_list)

diff_gene_vec <- paste0('sct_', diff_gene_vec)

test <- DoHeatmap(ilc.integrated,
                  features = diff_gene_vec,
                  group.by = "integrated_snn_res.0.4",
                  assay = "SCT",
                  slot = 'data')

###### Extended Data Fig 2d -Upset Plot for genes contained in the gene set list ######
# Copy the gene sets associated with each cluster
gene_sets_exp <- diff_gene_list

# Give cluster name to each geneset 
names(gene_sets_exp) <- c("Pre1","Pre2","Pre3", "Post1", "Post2","Post3")

# factor(gene_sets_exp, levels = c("Pre1","Pre2","Pre3", "Post1", "Post2","Post3"))

# Upset Plot
# https://cran.r-project.org/web/packages/UpSetR/vignettes/basic.usage.html

# Feature request to have features in each intersection set:
# https://github.com/hms-dbmi/UpSetR/issues/85

## Empty Set plot
# upset(fromList(gene_sets_exp), nsets = 6,
#       sets = c("Post3", "Post2","Post1","Pre3","Pre2","Pre1"),
#       order.by = c("freq"),
#       keep.order = TRUE,
#       empty.intersections = "on")

pdf(paste0(out.dir, "geneSet_RNA_clustering_upsetPlot_",exp.rep[i], ".pdf"),
    onefile=FALSE)

upset(fromList(gene_sets_exp), nsets = 6,
      sets = c("Post3", "Post2","Post1","Pre3","Pre2","Pre1"),
      order.by = c("freq"),
      keep.order = TRUE)

dev.off()


##### FIG. 2e - Enrichr Analysis #####
# Method:
# 1. Identify differential genes for each clsuter (padj < 0.05)
# 2. Run each gene set through Enrichr
# 3. Create Dot plot
# Note: Used Pyhton script to create text file for this analysis
# serody_project/bin/pre_post_ILC2_RNA_enrichr_analysis.ipynb

enrichr_path = '/proj/dllab/jfoster/serody_project/results/Fig2_mILC2_mILC2-GFP_scMultiome_analysis_RNA_Integration/enrichrData/hold/'

enrichr_data <- read.table(file = sprintf("%sILC2_enrichr_data.txt", enrichr_path),
           sep = '\t',
           quote = "",
           header = T)

enrichr_data$cluster <- as.factor(enrichr_data$cluster)

enrichr_data$cluster <- factor(enrichr_data$cluster, levels=c('clusterPost3',
                                                              'clusterPost2',
                                                              'clusterPost1'))

enrichr_plot <- ggplot(enrichr_data, aes(x=reorder(Term,-value,FUN=mean),
                         y=cluster,
                         size=value,
                         color=overlap )) + 
  geom_point() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90)) + 
  ylab("") + 
  xlab("")

enrichr_plot

ggsave(enrichr_plot, file = paste0(out.dir,"post_transplant_geneSet_Enrichr_MouseGeneAtlas_",exp.rep[i], ".pdf"),
       height = 7, width = 4 , device = "pdf")

##### FIG. 2h - Heatmap of differential chromatin accessibility for each cluster #####
# Note: Variable names have been carried over from the RNA abundance analysis.

###### Differential gene expression analysis of clusters ######

# Select clusters for a given resolution value
Idents(object = ilc.integrated) <- "integrated_snn_res.0.4"

# Store differential regions of chromatin accessibility for each cluster
diff_gene_list <- list()

# pre-transplant differential ATAC peaks relative to post-transplant cells
pre.transplant.clusters <- c(0,2,1)

for (pre.clust in pre.transplant.clusters) {
  
  # Identify differential ATAC peaks
  diff.genes <- FindMarkers(ilc.integrated,
                            ident.1 = pre.clust,
                            ident.2 = c(3,4,5),
                            min.pct = 0.25,
                            assay = "ATAC",
                            slot = "data")  
  
  # Select significant ATAC peaks with FC > 0
  top.clusterX.genes <- diff.genes %>%
    dplyr::filter(p_val_adj < 0.05 & avg_log2FC > 0) %>%
    arrange(desc(avg_log2FC))
  
  diff_gene_list <- append(diff_gene_list,list(row.names(top.clusterX.genes)))
  
  # # Write txt file of differential genes for each cluster
  # write.table(top.clusterX.genes, 
  #             paste0(out.dir, "ILC2_transplant_differential_genes_cluster_",pre.clust,".txt"),
  #             sep = '\t',
  #             row.names = T,
  #             quote = F)
  
} 


# post-transplant differential AAC peaks relative to all pre-transplant cells
post.transplant.clusters <- c(3,4,5)

for (post.clust in post.transplant.clusters) {
  
  # Identify differential genes
  diff.genes <- FindMarkers(ilc.integrated,
                            ident.1 = post.clust,
                            ident.2 = c(0,1,2),
                            min.pct = 0.25,
                            assay = "ATAC",
                            slot = "data")  
  
  # Select significant ATAC peaks with FC > 0
  top.clusterX.genes <- diff.genes %>%
    dplyr::filter(p_val_adj < 0.05 & avg_log2FC > 0) %>%
    arrange(desc(avg_log2FC))
  
  diff_gene_list <- append(diff_gene_list, list(row.names(top.clusterX.genes)))
  
  # Write txt file of differential genes for each cluster
  # write.table(top.clusterX.genes,paste0(out.dir, "ILC2_transplant_differential_genes_cluster_",post.clust,".txt"),
  #             sep = '\t',
  #             row.names = T,
  #             quote = F)
} 

###### Generate Union Set of differential genes  ######

# Clusters to cycle through
cluster.names <- levels(ilc.integrated)

cluster.names

#Need vector of cell IDs for each cluster 
cluster.cellIDs <- list()

#List of matricies of average matricies 
avg_gene_exp <- list()

# List of avg gene dataframes
gene_set_avg <- list()

for (gene_set in diff_gene_list) {
  
  for (cluster in cluster.names) {
    
    # Isolate the cellIDs for each cluster
    cell_IDs <- row.names(ilc.integrated@meta.data %>% dplyr::filter(integrated_snn_res.0.4 == cluster))
    
    # Normalized ATAC data
    rna_data <- ilc.integrated@assays[["ATAC"]]@data
    
    # Submatrix of normalized RNA signal
    sub_rna_data <- rna_data[gene_set, cell_IDs]
    
    # Vector of avg ATAC signal
    sub_rna_avg <- rowMeans(sub_rna_data)
    
    # Store avg ATAC signal (across cells) for each cluster
    avg_gene_exp <- append(avg_gene_exp, list(sub_rna_avg))
    
  }
  
  #Concat avgerage ATAC signal vectors into dataframe for each ATAC peak set
  gene_set_avg_df <-as.data.frame(do.call(cbind, avg_gene_exp))
  
  # Save avg. ATAC signal for each ATAC peak set for each cluster
  gene_set_avg <- append(gene_set_avg, list(gene_set_avg_df))
  
  # Reset the list holding the average score for each gene set
  avg_gene_exp <- list()
  
}

# Full Matrix of ATAC peaks 
diff_gene_clust_data <- do.call(rbind, gene_set_avg)

###### Complex Heatmap ######
# Mean Center the rows 
df <- t(scale(t(data.matrix(diff_gene_clust_data)), scale = FALSE))

#Save avg. chromatin accessibility dataFrame
# saveRDS(df, file = paste0(out.dir, "chromAccessibility_Heatmap_data_frame_",exp.rep[i], ".rds"))

# Read avg. chromatin accessibility dataFrame 
df <- readRDS(paste0(out.dir, "chromAccessibility_Heatmap_data_frame_rep2.rds"))

# Remove mean attribute from matrix - Not necessary 
attr(df, "scaled:center") <- NULL

cluster_names = c("Pre1","Pre2","Pre3",
                  "Post1","Post2","Post3")

row_split = rep(cluster_names, times = lengths(diff_gene_list))

#order row split with factor
row_split <- factor(row_split, levels=cluster_names)

pdf(paste0(out.dir, "pre_post_cluster_ATAC_heatmap_",exp.rep[i], ".pdf"),
    width=3,
    height=3.5)

Heatmap(df,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        row_names_gp = gpar(fontsize = 0),
        name = "Norm.\nGene Exp.",
        left_annotation = HeatmapAnnotation(clusters = row_split,
                                            col = list(clusters = c("Pre1" = "#005ba3",
                                                                    "Pre2" = "#a77f48",
                                                                    "Pre3" = "#63c7c5",
                                                                    "Post1" = "#238755",
                                                                    "Post2" = "#88b571",
                                                                    "Post3" = "#57813d")),
                                            which = 'row'),
        row_split = row_split,
        row_title = NULL,
        col = atac_color)


dev.off()

#### Fig2I-M - Putative TF regulators for pre- and post-transplant clusters ####
#Notes: 
# TF candidate analysis where comparison is only on the pre- or post-transplant
# groups.

# Select clustering resolution
resolution <- "integrated_snn_res.0.4"

# Set identity class to resolution 0.4 annotations
Idents(object = ilc.integrated) <- resolution

# Store cluster names to iterate through
clusters.names <- levels(ilc.integrated)

# Store number of motifs at gained regions of chromatin accessibility
num.motifs.at.gained <- list()

# Store number of putative regulators identified for each cluster 
num.putative.regulators <- list()

for (cluster in clusters.names) {
  
  ##### Putative regulators of pre-transplant clusters #####
  if (cluster == 0 | cluster == 1 | cluster == 2) {
   
    # Differential RNA analysis 
    markers_rna <- presto:::wilcoxauc.Seurat(X = ilc.integrated,
                                             group_by = resolution,
                                             assay = 'data',
                                             seurat_assay = 'SCT',
                                             groups_use = c(cluster, "3","4","5"))
    
    # Differential motif analysis 
    markers_motifs <- presto:::wilcoxauc.Seurat(X = ilc.integrated,
                                                group_by = resolution,
                                                assay = 'data',
                                                seurat_assay = 'chromvar',
                                                groups_use = c(cluster, "3","4","5"))
    
    # Store motif names
    motif.names <- markers_motifs$feature
    
    # Annotate differential RNA data with "RNA. "prefix - track RNA information 
    colnames(markers_rna) <- paste0("RNA.", colnames(markers_rna))
    
    # Annotate differential motif data with "motif. "prefix - track motif information
    colnames(markers_motifs) <- paste0("motif.", colnames(markers_motifs))
    
    # Add column with gene names
    markers_rna$gene <- markers_rna$RNA.feature
    
    # Add column with understandable motif name 
    markers_motifs$gene <- ConvertMotifID(ilc.integrated, id = motif.names, assay = "ATAC")
    
    # Function to identify candidate regulators 
    # TF with increased RNA and whose motif is present
    topTFs <- function(celltype, padj.cutoff = .05) {
      ctmarkers_rna <- dplyr::filter(
        markers_rna, RNA.group == celltype, RNA.padj < padj.cutoff, RNA.logFC > 0) %>% 
        arrange(-RNA.auc)
      ctmarkers_motif <- dplyr::filter(
        markers_motifs, motif.group == celltype, motif.padj < padj.cutoff, motif.logFC > 0) %>% 
        arrange(-motif.auc)
      
      ctmarkers_motif$gene <- tolower(ctmarkers_motif$gene)
      
      ctmarkers_motif$gene <- str_to_title(ctmarkers_motif$gene)
      
      top_tfs <- inner_join(
        x = ctmarkers_rna[, c(2, 11, 6, 7)], 
        y = ctmarkers_motif[, c(2, 1, 11, 6, 7)], by = "gene"
      )
      top_tfs$avg_auc <- (top_tfs$RNA.auc + top_tfs$motif.auc) / 2
      top_tfs <- arrange(top_tfs, -avg_auc)
      return(top_tfs)
    }
    
    # Top 10 putative regulator for a given cluster
    clusterX.regulators <- head(topTFs(cluster), n = 10)
    
    # Store number of putative regulator for a given cluster
    num.putative.regulators <- append(num.putative.regulators, list(nrow(topTFs(cluster))))
    
    # Bar plot of top 10 putative regulator for a given cluster
    clusterX.regulators %>% 
      ggplot(aes(x = avg_auc,y = gene)) +
      geom_bar(stat = "identity")+
      scale_fill_gradient()+
      aes(y = fct_reorder(gene,avg_auc)) +
      theme_classic() +
      xlab("Average AUC") + 
      ylab("Candiate Regulators") +
      ggtitle(paste0("Cluster ",cluster)) + 
      scale_x_continuous(limits = c(0,1), breaks= c(0, 0.5, 1), labels=c(0, 0.5, 1)) +
      geom_text(aes(label= round(avg_auc, digits = 2)), vjust=1.25, hjust=1.25, color="white", size=3.5)
    
    output.filename = paste0(out.dir, "cluster",cluster,"_v_postClusters_candidate_TF_",exp.rep[i],".pdf")
    
    ggsave(file = output.filename,height = 4.6, width = 3, device = "pdf")
  }
  
  ##### Putative regulators of post-transplant clusters #####
  else if (cluster == 3 | cluster == 4 | cluster == 5) {
   
    # Differential RNA analysis 
    markers_rna <- presto:::wilcoxauc.Seurat(X = ilc.integrated,
                                             group_by = resolution,
                                             assay = 'data',
                                             seurat_assay = 'SCT',
                                             groups_use = c(cluster, "0","1","2"))
    
    # Differential motif analysis 
    markers_motifs <- presto:::wilcoxauc.Seurat(X = ilc.integrated,
                                                group_by = resolution,
                                                assay = 'data',
                                                seurat_assay = 'chromvar',
                                                groups_use = c(cluster, "0","1","2"))
    
    # Store motif names
    motif.names <- markers_motifs$feature
    
    # Annotate differential RNA data with "RNA. "prefix - track RNA information 
    colnames(markers_rna) <- paste0("RNA.", colnames(markers_rna))
    
    # Annotate differential motif data with "motif. "prefix - track motif information
    colnames(markers_motifs) <- paste0("motif.", colnames(markers_motifs))
    
    # Add column with gene names
    markers_rna$gene <- markers_rna$RNA.feature
    
    # Add column with understandable motif name 
    markers_motifs$gene <- ConvertMotifID(ilc.integrated, id = motif.names, assay = "ATAC")
    
    # Function to identify candidate regulators 
    # TF with increased RNA and whose motif is present
    topTFs <- function(celltype, padj.cutoff = .05) {
      ctmarkers_rna <- dplyr::filter(
        markers_rna, RNA.group == celltype, RNA.padj < padj.cutoff, RNA.logFC > 0) %>% 
        arrange(-RNA.auc)
      ctmarkers_motif <- dplyr::filter(
        markers_motifs, motif.group == celltype, motif.padj < padj.cutoff, motif.logFC > 0) %>% 
        arrange(-motif.auc)
      
      ctmarkers_motif$gene <- tolower(ctmarkers_motif$gene)
      
      ctmarkers_motif$gene <- str_to_title(ctmarkers_motif$gene)
      
      top_tfs <- inner_join(
        x = ctmarkers_rna[, c(2, 11, 6, 7)], 
        y = ctmarkers_motif[, c(2, 1, 11, 6, 7)], by = "gene"
      )
      top_tfs$avg_auc <- (top_tfs$RNA.auc + top_tfs$motif.auc) / 2
      top_tfs <- arrange(top_tfs, -avg_auc)
      return(top_tfs)
    }
    
    # Top 10 putative regulator for a given cluster
    clusterX.regulators <- head(topTFs(cluster),n = 10)
    
    # Store number of putative regulator for a given cluster
    num.putative.regulators <- append(num.putative.regulators, list(nrow(topTFs(cluster))))
    
    # Bar plot of top 10 putative regulator for a given cluster
    clusterX.regulators %>% 
      ggplot(aes(x = avg_auc,y = gene)) +
      geom_bar(stat = "identity")+
      scale_fill_gradient()+
      aes(y = fct_reorder(gene,avg_auc)) +
      theme_classic() +
      xlab("Average AUC") + 
      ylab("Candiate Regulators") +
      ggtitle(paste0("Cluster ",cluster)) + 
      scale_x_continuous(limits = c(0,1), breaks= c(0, 0.5, 1), labels=c(0, 0.5, 1)) +
      geom_text(aes(label= round(avg_auc, digits = 2)), vjust=1.25, hjust=1.25, color="white", size=3.5)
      
    
    output.filename = paste0(out.dir, "cluster",cluster,"_v_postClusters_candidate_TF_",exp.rep[i],".pdf")
    
    ggsave(file = output.filename,height = 4.6, width = 3, device = "pdf")  
  }
  
  ##### Number of differential motifs at gained regions ####
  
  # Number of differential motifs at gained regions of chromatin accessibility
  num.motifs.cluster <- markers_motifs %>% 
    dplyr::filter(motif.group == cluster & motif.logFC > 0 & -log10(motif.padj) > -log10(0.05)) %>%
    arrange(desc(motif.logFC)) %>% nrow()
  
  num.motifs.at.gained <- append(num.motifs.at.gained, list(num.motifs.cluster))
  
  ##### Volcano plot of differential motifs #####
  
  # Direction of Log2Fold change is relative to a particular cluster 
  # Here I've selected either the pre- or post- cluster that is be tested
  # against the entirety of post- or pre- samples respectively
  motif.data <- markers_motifs %>% 
    dplyr::filter(motif.group == cluster) %>% 
    arrange(desc(motif.logFC))
  
  # Add column that specifies the color of significant or non-significant 
  # motifs
  motif.volcano.data <- motif.data %>% 
    mutate(sig.color  = ifelse(-log10(motif.padj) > -log10(0.05), '#1F9FE1', '#E0601E')) %>%
    arrange(sig.color)
  
  # Volcano Plot - highlighting specific ILC1 and ILC2 genes
  ggplot(motif.volcano.data, aes(x=motif.logFC, y=-log10(motif.padj),  col = sig.color, label = gene)) + 
    geom_point() +
    scale_color_identity() +
    geom_text_repel(box.padding = .5,
                    size=3,
                    segment.size = 0.25,
                    force_pull=1,
                    max.overlaps=50,
                    min.segment.length= 0.01,
                    data=subset(motif.ma.data, (gene == "GATA3" | gene == "TBX21" | gene == "BATF" |
                                                  gene == "MAF" | gene == "FOSB" | gene == "BACH1" |
                                                  gene == "EOMES" | gene == "RBPJL" | gene == "FOXO3" |
                                                  gene == "TCF3" | gene == "SMAD2" | gene == "BATF" |
                                                  gene == "NFATC4" | gene == "NR4A1" | gene == "CREB1" |
                                                  gene == "MECOM" | gene == "IRF8" | gene == "NR4A2" |
                                                  gene == "POU2F2" | gene == "SNAI3")),
                    colour = 'black') +
    theme_classic() + 
    theme(legend.position="none") +
    ylab("-log10(p-value)") +
    xlab(paste0("Motif Enrichment \n log2(Cluster",cluster, "/ all)"))
  
  # Save Volcano plot
  ggsave(file = paste0(out.dir, "pre_post_diff_motif_volcano_cluster_",cluster,"_",exp.rep[i],".pdf"),
         height = 4.5, width = 4.5, device = "pdf")
  
}

