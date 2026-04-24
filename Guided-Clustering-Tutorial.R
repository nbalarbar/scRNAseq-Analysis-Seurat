library(languagserver)
library(dplyr)
library(Seurat)
library(patchwork)
# ===========================================================================
# DATA SETUP AND SEURAT OBJECT CREATION
# ===========================================================================

# Load Dataset
pbmcData <- Read10X(data.dir = "data/filtered_gene_bc_matrices/hg19/")
# Initialize Seurat object with raw data
# pbmc objects contain data, analysis for single cell data sets
pbmc <- CreateSeuratObject(counts = pbmcData, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc

# Output:
# An object of class Seurat 
# 13714 features across 2700 samples within 1 assay 
# Active assay: RNA (13714 features, 0 variable features)
#  1 layer present: counts


# Look at the first few genes in the first 30 cells
pbmcData[c("CD3D", "TCL1A", "MS4A1"), 1:30]
# Columns: Discrete samples
# Rows: Given Gene
# Value: How many molecules detected (. == 0)

#3 x 30 sparse Matrix of class "dgCMatrix"                       
# CD3D  4 . 10 . . 1 2 3 1 . . 2 7 1 . . 1 3 . 2  3 . . . . . 3 4 1 5
# TCL1A . .  . . . . . . 1 . . . . . . . . . . .  . 1 . . . . . . . .
# MS4A1 . 6  . . . . . . 1 1 1 . . . . . . . . . 36 1 2 . . 2 . . . .


# ===========================================================================
# PRE-PROCESSING WORKFLOW
# ===========================================================================

# QC AND SELECTING VIABLE CELLS

# The [[operator addes columns to meta data]], Stash QC stats here

pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, patter = "^MT-") # Track Mitochondria distress

# Check QC metrics
head(pbmc@meta.data, 5)

# Output:
# > head(pbmc@meta.data, 5)
#                  orig.ident nCount_RNA nFeature_RNA percent.mt
# AAACATACAACCAC-1     pbmc3k       2419          779  3.0177759
# AAACATTGAGCTAC-1     pbmc3k       4903         1352  3.7935958
# AAACATTGATCAGC-1     pbmc3k       3147         1129  0.8897363
# AAACCGTGCTTCCG-1     pbmc3k       2639          960  1.7430845
# AAACCGTGTATGCG-1     pbmc3k        980          521  1.2244898

# Note the new "percent.mt" column

# We'll only keep cells that have <5% mitochondrial counts
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# FeatureScatter is used to visualize feature-feature relationships:
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA <2500 & percent.mt < 5)

# NORMALIZE DATA
# By Default, we use LogNormalize, then multiply by a scale factor of 10,000

# IMPORTANT: Normalized data is stored in pbmc[["RNA"]]$data

pbmc <- NormalizeData(pbmc)
# Default Parameters:
# pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)


# FEATURE SELECTION
# Identify highly variable features, by default 2000
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Take top 10 most variable genes
top10 <- head(VariableFeatures(pbmc), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# SCALE DATA
# Needed prior to PCA
# Results are stored in pbmc[["RNA"]]$scale.data

allGenes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = allGenes)
