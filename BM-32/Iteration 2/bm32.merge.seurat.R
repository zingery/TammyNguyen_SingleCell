library(tidyverse) 
library(data.table)
library(Seurat)


##------ Create seurat object ------

# Read 10x output
bmy32.10x<- Read10X(data.dir = "~/Dropbox (UMass Medical School)/Zinger_at_UMass/Silvia Coreva/Tammy/singleCell/BM32_SingleCell/Iteration1_CellRanger/BMY-32/")
bmr32.10x<- Read10X(data.dir = "~/Dropbox (UMass Medical School)/Zinger_at_UMass/Silvia Coreva/Tammy/singleCell/BM32_SingleCell/Iteration1_CellRanger/BMR-32/")
bms32.10x<- Read10X(data.dir = "~/Dropbox (UMass Medical School)/Zinger_at_UMass/Silvia Coreva/Tammy/singleCell/BM32_SingleCell/Iteration1_CellRanger/BMS-32/")

bmy32.seurat <- CreateSeuratObject(counts = bmy32.10x, project = "BMY-32" ) 
bmr32.seurat <- CreateSeuratObject(counts = bmr32.10x, project = "BMR-32" ) 
bms32.seurat <- CreateSeuratObject(counts = bms32.10x, project = "BMS-32" ) 


# Merge datasets into one single seurat object
bm32.merge.seurat <- merge(bmy32.seurat, c(bmr32.seurat, bms32.seurat), 
                 add.cell.ids = c("BMY-32", "BMR-32", "BMS-32"))

rm(bmy32.10x,bmr32.10x,bms32.10x)
rm(bmy32.seurat,bmr32.seurat,bms32.seurat)
gc()

# Number of cells in each sample in this seurat object
table(bm32.merge.seurat$orig.ident)


## Normalization
bm32.merge.seurat <- SCTransform(bm32.merge.seurat, assay = 'RNA')

## PCA
bm32.merge.seurat <- RunPCA(object = bm32.merge.seurat ) 

mat <- Seurat::GetAssayData(bm32.merge.seurat, assay = "SCT", slot = "scale.data")
pca <- undiff.seurat[["pca"]]
total_variance <- sum(matrixStats::rowVars(mat))   # total variance
eigValues = (pca@stdev)^2  # EigenValues
varExplained = eigValues / total_variance
plot(c(1:50),varExplained) # elbow plot

## UMAP
bm32.merge.seurat <- FindNeighbors(object = bm32.merge.seurat)
bm32.merge.seurat <- RunUMAP(object = bm32.merge.seurat, dims = 1:20) # dims approximately determined by elbow plot

## Clustering (higher resolution = more number of clusters)
bm32.merge.seurat <- FindClusters(object = bm32.merge.seurat, resolution = 0.2)

# How many cells are in each cluster
table(Idents(bm32.merge.seurat))

saveRDS(bm32.merge.seurat, file="~/Dropbox (UMass Medical School)/Zinger_at_UMass/Silvia Coreva/Tammy/singleCell/BM32_SingleCell/Iteration2_SeuratNoFilter/bm32.merge.seurat.RDS")


##-----------------------------------------------------------------------------------------------



bm32.merge.seurat <- readRDS("~/Dropbox (UMass Medical School)/Zinger_at_UMass/Silvia Coreva/Tammy/singleCell/BM32_SingleCell/Iteration2_SeuratNoFilter/bm32.merge.seurat.RDS")

##----- Plot -----

# UMAP, color by cluster
DimPlot(bm32.merge.seurat, reduction = "umap", label=T ) # with cluster label

# UMAP, color by sample ID
DimPlot(bm32.merge.seurat, reduction = "umap", group.by = "orig.ident", cols=c("#d62728","#1f77b4","#e7ba52"))  

# Plot gene expression on UMAP projection
FeaturePlot(bm32.merge.seurat, features = c("WNT5A"), ncol = 1, 
            cols=c("#eeeeee","#d62728"),  
            pt.size = 0.5, order=T) & NoAxes()

 
##----- QC -----
# https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_01_qc.html
bm32.merge.seurat <- PercentageFeatureSet(bm32.merge.seurat, "^MT-", col.name = "percent_mito")
bm32.merge.seurat <- PercentageFeatureSet(bm32.merge.seurat, "^RP[SL]", col.name = "percent_ribo")
bm32.merge.seurat <- PercentageFeatureSet(bm32.merge.seurat, "^HB[^(P)]", col.name = "percent_hemoglobin")
                      # hemoglobin genes - includes all genes starting with HB except HBP.


# Violin plot by cluster
VlnPlot(bm32.merge.seurat,  
        features = c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "percent_hemoglobin"),  
        pt.size = 0, ncol = 5) + 
  NoLegend()

# Violin plot by sample
VlnPlot(bm32.merge.seurat, group.by = "orig.ident", 
        features = c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "percent_hemoglobin"), 
        cols=c("#dedede","#1f77b4","#dedede") ,
        pt.size = 0, ncol = 5) + 
  NoLegend()


##----- Differential Expression ------
# Cluster markers 
# (ident.1 cluster vs all others, genes must be detected in at least 25% of cells in one of the comparison group)
# (log2 fold change >1 or <-1 results are returned)
marker.clust9 <- FindMarkers(bm32.merge.seurat, ident.1 = "9", min.pct = 0.25, logfc.threshold = 1)
marker.clust8 <- FindMarkers(bm32.merge.seurat, ident.1 = "8", min.pct = 0.25, logfc.threshold = 1) 
marker.clust7 <- FindMarkers(bm32.merge.seurat, ident.1 = "7", min.pct = 0.25, logfc.threshold = 1) 
 
# Comparing two clusters 
marker.clust0v5 <- FindMarkers(bm32.merge.seurat, ident.1 = "0",ident.2 = "5", min.pct = 0.25, logfc.threshold = 1)
 