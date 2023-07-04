# Calculate the gene-level RNA abundance estimates from transcript-level
# salmon files. 


library(tximport)

# load RefSeq Gene Data Table
tx2gene <- read.table("./referenceData/RefSeq_mm10_030617.txt", header=F, sep="\t")

# allow command line arguments
args = commandArgs(TRUE)

# first argument is table of sample names with correpsonding paths to salmon .sf files
samples=read.table("./Figure_1/bruce_sf_files_formatted.txt",header=F,row.names=1)

# makes a file object of all paths to sf files
files <- file.path(samples$V2)

# assigns sample name to each file pathe of file object
names(files) <- rownames(samples)

# convert transcript level summaries to gene level summaries
txi <- tximport(files,
                type = "salmon",
                tx2gene = tx2gene,
                importer = read.delim,
                txOut=F)

# Output gene-level data
write.table(txi,
            file="./Figure_1/ILC2_Bruce_geneLevel_RNA_abundance.txt",
            col.names=NA,
            row.names=T,
            sep="\t",
            quote=F)

# Save SessionInfo 
writeLines(capture.output(sessionInfo()), "./sessionInfo/transcritToGeneConversion_mm10_sessionInfo.txt")
