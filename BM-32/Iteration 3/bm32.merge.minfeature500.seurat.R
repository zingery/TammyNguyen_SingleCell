library(tidyverse) 
library(data.table)
library(Seurat)


# Create seurat object
bmy32.10x<- Read10X(data.dir = "~/Dropbox (UMass Medical School)/Zinger_at_UMass/Silvia Coreva/Tammy/singleCell/BM32_SingleCell/Iteration1_CellRanger/BMY-32/")
bmr32.10x<- Read10X(data.dir = "~/Dropbox (UMass Medical School)/Zinger_at_UMass/Silvia Coreva/Tammy/singleCell/BM32_SingleCell/Iteration1_CellRanger/BMR-32/")
bms32.10x<- Read10X(data.dir = "~/Dropbox (UMass Medical School)/Zinger_at_UMass/Silvia Coreva/Tammy/singleCell/BM32_SingleCell/Iteration1_CellRanger/BMS-32/")

# create object with minimum feature requirement
bmy32.seurat <- CreateSeuratObject(counts = bmy32.10x, project = "BMY-32", min.features = 500, min.cells = 10 ) 
bmr32.seurat <- CreateSeuratObject(counts = bmr32.10x, project = "BMR-32", min.features = 500, min.cells = 10 )  
bms32.seurat <- CreateSeuratObject(counts = bms32.10x, project = "BMS-32", min.features = 500, min.cells = 10 )  

# Merge datasets into one single seurat object
bm32.minfeature500.seurat <- merge(bmy32.seurat, c(bmr32.seurat, bms32.seurat), 
                 add.cell.ids = c("BMY-32", "BMR-32", "BMS-32"))

table(bm32.minfeature500.seurat$orig.ident)
# BMR-32 BMS-32 BMY-32 
# 16183   9095   6981 

dim(bm32.minfeature500.seurat)
# [1] 21998 32259

rm(bmy32.10x,bmr32.10x,bms32.10x)
rm(bmy32.seurat,bmr32.seurat,bms32.seurat)
gc()

## Normalization round 1
bm32.minfeature500.seurat <- SCTransform(bm32.minfeature500.seurat, assay = 'RNA')

## PCA round 1
bm32.minfeature500.seurat <- RunPCA(object = bm32.minfeature500.seurat ) 


mat <- Seurat::GetAssayData(bm32.minfeature500.seurat, assay = "SCT", slot = "scale.data")
pca <- bm32.minfeature500.seurat[["pca"]]
total_variance <- sum(matrixStats::rowVars(mat))   # Get the total variance:
eigValues = (pca@stdev)^2  ## EigenValues
varExplained = eigValues / total_variance

plot(c(1:50),varExplained)

## UMAP   round 1
bm32.minfeature500.seurat <- FindNeighbors(object = bm32.minfeature500.seurat)
bm32.minfeature500.seurat <- FindClusters(object = bm32.minfeature500.seurat, resolution = 0.2)
bm32.minfeature500.seurat <- RunUMAP(object = bm32.minfeature500.seurat, dims = 1:20)

##### filter and re-normalize/cluster/etc
# Remove cluster 6 and cluster 8 (poor nuclear lysis, high mito)
bm32.minfeature500.seurat <- subset(bm32.minfeature500.seurat, idents = c('6','8'), invert=TRUE)

# How many cells are in each sample
table(bm32.minfeature500.seurat$orig.ident)
# BMR-32 BMS-32 BMY-32 
# 15702   8951   5945 

# How many cells are in each cluster
table(Idents(bm32.minfeature500.seurat))

## Normalization after filter
bm32.minfeature500.seurat <- SCTransform(bm32.minfeature500.seurat, assay = 'RNA')

## PCA after filter
bm32.minfeature500.seurat <- RunPCA(object = bm32.minfeature500.seurat ) 


mat <- Seurat::GetAssayData(bm32.minfeature500.seurat, assay = "SCT", slot = "scale.data")
pca <- bm32.minfeature500.seurat[["pca"]]
total_variance <- sum(matrixStats::rowVars(mat))   # Get the total variance:
eigValues = (pca@stdev)^2  ## EigenValues
varExplained = eigValues / total_variance

plot(c(1:50),varExplained)

## UMAP
bm32.minfeature500.seurat <- FindNeighbors(object = bm32.minfeature500.seurat)
bm32.minfeature500.seurat <- FindClusters(object = bm32.minfeature500.seurat, resolution = 0.5)
bm32.minfeature500.seurat <- RunUMAP(object = bm32.minfeature500.seurat, dims = 1:20)

saveRDS(bm32.minfeature500.seurat, file="~/Dropbox (UMass Medical School)/Zinger_at_UMass/Silvia Coreva/Tammy/singleCell/BM32_SingleCell/Iteration3_SeuratMinFeature/bm32.minfeature500.seurat.RDS")




## Plot PCA, UMAP
DimPlot(bm32.minfeature500.seurat, reduction = "pca" )
DimPlot(bm32.minfeature500.seurat, reduction = "pca" , group.by = "orig.ident", cols=c("#d62728","#1f77b4","#e7ba52") ) # color cells by sample

DimPlot(bm32.minfeature500.seurat, reduction = "umap", label=T ) # label the clusters
DimPlot(bm32.minfeature500.seurat, reduction = "umap", group.by = "orig.ident", cols=c("#d62728","#1f77b4","#e7ba52") ) # color cells by sample



## Plot gene expression on UMAP projection
FeaturePlot(bm32.minfeature500.seurat, features = c("LEPR"), ncol = 1, 
            cols=c("#eeeeee","#d62728"),  
            pt.size = 0.5, order=T) & NoAxes()

 
## Export UMAP coordinates for Loupe browser
bm32.minfeature500.seurat.umap.coordinates <- as.data.frame(Embeddings(object = bm32.minfeature500.seurat[["umap"]])) %>%
  tibble::rownames_to_column(var = "sample_barcode") %>%
  separate(sample_barcode, into=c("sample","barcode"), sep="_") %>%
  separate(barcode, into=c("barcode1","sample1"), sep="-") %>%
  mutate(sample1 = case_when(
    sample == "BMY-32" ~ "1",
    sample == "BMR-32" ~ "2",
    sample == "BMS-32" ~ "3"
  )) %>%
  mutate(barcode = paste0(barcode1,"-",sample1)) %>%
  select(barcode, UMAP_1,UMAP_2)

write.csv(bm32.minfeature500.seurat.umap.coordinates, file="~/Documents/GitHub/TammyNguyen_SingleCell/BM-32/Iteration 3/bm32.minfeature500.seurat.umap.csv", row.names=F, quote=F)
  
 
#bm32.minfeature500.seurat <- FindClusters(object = bm32.minfeature500.seurat, resolution = 0.425)  

# resolution
# 0.1: 5 clusters
# 0.12: 6 clusters
# 0.15: 7 clusters
# 0.2: 8 clusters
# 0.25: 9 clusters
# 0.3: 10 clusters
# 0.35: 11 clusters
# 0.4: 12 clusters
# 0.45: 15 clusters
# 0.5: 16 clusters

bm32.minfeature500.seurat.clusters <- bm32.minfeature500.seurat@meta.data %>% 
  as.data.frame() %>%
  tibble::rownames_to_column(var = "sample_barcode") %>%
  separate(sample_barcode, into=c("sample","barcode"), sep="_") %>%
  separate(barcode, into=c("barcode1","sample1"), sep="-") %>%
  mutate(sample1 = case_when(
    sample == "BMY-32" ~ "1",
    sample == "BMR-32" ~ "2",
    sample == "BMS-32" ~ "3"
  )) %>%
  mutate(barcode = paste0(barcode1,"-",sample1)) %>%
  select(barcode,`SCT_snn_res.0.1`, `SCT_snn_res.0.12`,`SCT_snn_res.0.15`,
         `SCT_snn_res.0.2`,`SCT_snn_res.0.25`,`SCT_snn_res.0.3`,`SCT_snn_res.0.35`,`SCT_snn_res.0.4`)


colnames(bm32.minfeature500.seurat.clusters) <-
  c("barcode","Seurat_min.feature500_5clusters","Seurat_min.feature500_6clusters",
    "Seurat_min.feature500_7clusters","Seurat_min.feature500_8clusters",
    "Seurat_min.feature500_9clusters","Seurat_min.feature500_10clusters",
    "Seurat_min.feature500_11clusters","Seurat_min.feature500_12clusters")

write.csv(bm32.minfeature500.seurat.clusters, file="~/Documents/GitHub/TammyNguyen_SingleCell/BM-32/Iteration 3/bm32.merge.minfeature500.seurat.clusters.csv", row.names=F, quote=F)

 

##----- QC -----
# https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_01_qc.html
bm32.minfeature500.seurat <- PercentageFeatureSet(bm32.minfeature500.seurat, "^MT-", col.name = "percent_mito")
bm32.minfeature500.seurat <- PercentageFeatureSet(bm32.minfeature500.seurat, "^RP[SL]", col.name = "percent_ribo")
bm32.minfeature500.seurat <- PercentageFeatureSet(bm32.minfeature500.seurat, "^HB[^(P)]", col.name = "percent_hemoglobin")
# hemoglobin genes - includes all genes starting with HB except HBP.


# Violin plot by cluster
VlnPlot(bm32.minfeature500.seurat,  
        features = c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "percent_hemoglobin"),  
        pt.size = 0, ncol = 5) + 
  NoLegend()

# Violin plot by sample
VlnPlot(bm32.minfeature500.seurat, group.by = "orig.ident", 
        features = c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "percent_hemoglobin"),  
        pt.size = 0, ncol = 5) + 
  NoLegend()


##----- Differential Expression ------
# Cluster markers 
# (ident.1 cluster vs all others, genes must be detected in at least 25% of cells in one of the comparison group)
# (log2 fold change >1 or <-1 results are returned)
marker.clust9 <- FindMarkers(bm32.minfeature500.seurat, ident.1 = "9", min.pct = 0.25, logfc.threshold = 1)
marker.clust8 <- FindMarkers(bm32.minfeature500.seurat, ident.1 = "8", min.pct = 0.25, logfc.threshold = 1) 
marker.clust7 <- FindMarkers(bm32.minfeature500.seurat, ident.1 = "7", min.pct = 0.25, logfc.threshold = 1) 

# Comparing two clusters 
marker.clust0v5 <- FindMarkers(bm32.minfeature500.seurat, ident.1 = "0",ident.2 = "5", min.pct = 0.25, logfc.threshold = 1)
