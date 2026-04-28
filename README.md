# **Cardiac scRNA-seq: Healthy vs. Myocardial Infarction**
**Single-Cell Transcriptomic Analysis of Healthy Cardiac and Acute Myocardial Infarction (AMI) cells**

**Status: Work in Progress**

**Tools:** R, Seuratv5, SQL

## **Overview**
This project investigates transcriptomic differences between healthy and acutely infarcted cardiac tissue at single-cell resolution. Using publicly available scRNA-seq data (GSE180678), I identify distinct cell populations, annotate cell types via marker gene expression, and quantify shifts in cell type composition between disease conditions.

This work serves as a personal bridge between my experience in clinical cardiovascualar research and my MSc in Bioinformatics and Computational Biology at The University of Texas at Dallas. Applying computational methods to single cell genomics to answer questions relating to acute cardiovascular events.

 This project aims to build directly on my prior work contributing to the DILWALE registry, the largest South Asian American EHR-derived cardiovascular dataset in Texas, moving from population-level clinical data toward molecular-level transcriptomic analysis.

**Central Question:** With respect to gene expression, what distinct cell populations exist in cardiac tissue, and how does acute myocardial infarction alter gene expression?

## **Dataset**
| Header 1 | Header 2 |
| -------- | -------- | 
| Source  | NCBI Gene Expression Omnibus|
| Acession ID  | GSE180678  |
| Organism  | Homosapien|
| Tissue | Cardiac  |
| Cells(Pre-QC)  | 1,290|
| Conditions  | Healthy vs AMI cells|

## **Methods/Pipeline**
1) Raw Count Matrix
2) Quality Control
    - (mitochondrial %, nFeature, nCount filtering)
3) SCTransform
4) PCA
    - Elbow Plot + PC Selection
5) Clustering
    - KNN
6) UMAP Cluster Visualization
7) Cell Type Annotation
8) Cell Type Comparision (Healthy vs. AMI)
9) Metadata Export
    - SQLite Database
    - SQL Queries
## **SQL Integration**
After annotation, all per-cell metadata is exported from Seurat and stored in a local SQLite database. This enables structured querying of biological results without re-running the full Seurat pipeline.

_In progress.._

## **Repository Structure**
_In Progress..._
## **Key Findings**
_In progress..._
## **Author**
Noah Balarbar MSc. Computational Biology and Bioinformatics, University of Texas at Dallas (Expected 2026) 

[LinkedIn](www.linkedin.com/in/noah-balarbar)

noah.balarbar@gmail.com