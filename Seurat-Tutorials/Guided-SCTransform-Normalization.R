library(Seurat)
library(ggplot2)
library(sctransform)

pbmcData <- Read10X("data/filtered_gene_bc_matrices/hg19")
pbmc <- CreateSeuratObject(counts = pbmcData)
