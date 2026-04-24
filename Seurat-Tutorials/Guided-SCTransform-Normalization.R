library(Seurat)
library(ggplot2)
library(sctransform)
# ===========================================================================
# DATA SETUP AND SEURAT OBJECT CREATION
# ===========================================================================
pbmcData <- Read10X("data/filtered_gene_bc_matrices/hg19")
pbmc <- CreateSeuratObject(counts = pbmcData)

# ===========================================================================
# APPLY SCTRANSFORM NORMALIZATION
# ===========================================================================
# SCTransform completes all 3: Normalization, Scaling, and Finding Variable Features all in a single command
# Replaces these steps in the general workflow
# In addition, SCTransform is also able to remove confounding sources of variation

# NOTE: This will run faster when "glmGamPoi" library is also installed

# Mitochondrial percentage
pbmc <- PercentageFeatureSet(pbmc, pattern = "^MT-", col.name = "percent.mt")

# Perform sctransform
pbmc <- SCTransform(pbmc, vars.to.regress = "percent.mt", verbose = FALSE)

# ===========================================================================
# PCA AND UMAP
# ===========================================================================
pbmc <- RunPCA(pbmc, verbose = FALSE)
pbmc <- RunUMAP(pbmc, dims = 1:30, verbose = FALSE)

# The following runs on PCA not UMAP
pbmc <- FindNeighbors(pbmc, dims = 1:30, verbose = FALSE)
pbmc <- FindClusters(pbmc, verbose = FALSE)

# UMAP Visualization
DimPlot(pbmc, label = TRUE)
