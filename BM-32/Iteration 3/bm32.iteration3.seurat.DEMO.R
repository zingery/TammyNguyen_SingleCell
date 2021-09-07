library(tidyverse) 
library(data.table)
library(Seurat)

 

##----- Read the pre-calculated Seurat object -----
# remember to change the file path to the path on your computer
bm32.minfeature2000.seurat <- readRDS("Iteration3_SeuratMinFeature/bm32.minfeature2000.seurat.RDS")


##----- Plot the clusters/samples in PCA, UMAP projections -----

# PCA, color by cluster 
DimPlot(bm32.minfeature2000.seurat, reduction = "pca", label=T )

# PCA, color by sample (original identifier)
DimPlot(bm32.minfeature2000.seurat, reduction = "pca" , group.by = "orig.ident", cols=c("#d62728","#1f77b4","#e7ba52") )

# UMAP, color by cluster 
DimPlot(bm32.minfeature2000.seurat, reduction = "umap", label=T )

# UMAP, color by sample (original identifier) 
DimPlot(bm32.minfeature2000.seurat, reduction = "umap", group.by = "orig.ident", cols=c("#d62728","#1f77b4","#e7ba52") )

# Plot a gene's expression on UMAP projection
FeaturePlot(bm32.minfeature2000.seurat, features = c("CEBPA"), order=T) 

# Plot a gene's expression on violin plot, separate by cluster
VlnPlot(bm32.minfeature2000.seurat,  features = c("CEBPA") )  



##----- Adjust Clustering Resolution -----
# adjust the "resolution" parameter for different number of clusters
# higher resolution = more number of clusters

bm32.minfeature2000.seurat <- FindClusters(object = bm32.minfeature2000.seurat, resolution = 0.11)
 
# bm32.minfeature2000.seurat resolutions
# 0.1: 5 clusters #
# 0.1075: 6 clusters #
# 0.11: 7 clusters #
# 0.2: 8 clusters #
# 0.36: 9 clusters # 
# 0.38: 10 clusters   
# 0.4: 12 clusters # 
# 0.5: 14 clusters 
 
 

##----- QC plots -----
# https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_01_qc.html

# Calculate the percentage of mitochondrial/ribosomal genes
bm32.minfeature2000.seurat <- PercentageFeatureSet(bm32.minfeature2000.seurat, "^MT-", col.name = "percent_mito")
bm32.minfeature2000.seurat <- PercentageFeatureSet(bm32.minfeature2000.seurat, "^RP[SL]", col.name = "percent_ribo")
 
# Violin plot by cluster
VlnPlot(bm32.minfeature2000.seurat,  
        features = c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo" ),  
        pt.size = 0, ncol = 4)  

 

##----- Differential Expression ------

# Identify cluster markers 
# (one cluster vs all others, genes must be detected in at least 25% of cells in one of the comparison group)
# (log2 fold change >1 or <-1 results are returned)
marker.clust9 <- FindMarkers(bm32.minfeature2000.seurat, ident.1 = "9", min.pct = 0.25, logfc.threshold = 1)
marker.clust8 <- FindMarkers(bm32.minfeature2000.seurat, ident.1 = "8", min.pct = 0.25, logfc.threshold = 1) 
marker.clust7 <- FindMarkers(bm32.minfeature2000.seurat, ident.1 = "7", min.pct = 0.25, logfc.threshold = 1) 

# Comparing two clusters 
marker.clust0v5 <- FindMarkers(bm32.minfeature2000.seurat, ident.1 = "0",ident.2 = "5", min.pct = 0.25, logfc.threshold = 1)
