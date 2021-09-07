library(tidyverse) 
library(data.table)
library(Seurat)


# Create seurat object
bmy32.10x<- Read10X(data.dir = "~/Dropbox (UMass Medical School)/Zinger_at_UMass/Silvia Coreva/Tammy/singleCell/BM32_SingleCell/Iteration1_CellRanger/BMY-32/")
bmr32.10x<- Read10X(data.dir = "~/Dropbox (UMass Medical School)/Zinger_at_UMass/Silvia Coreva/Tammy/singleCell/BM32_SingleCell/Iteration1_CellRanger/BMR-32/")
bms32.10x<- Read10X(data.dir = "~/Dropbox (UMass Medical School)/Zinger_at_UMass/Silvia Coreva/Tammy/singleCell/BM32_SingleCell/Iteration1_CellRanger/BMS-32/")

# create object with minimum feature requirement
bmy32.seurat <- CreateSeuratObject(counts = bmy32.10x, project = "BMY-32", min.features = 2000, min.cells = 10 ) 
bmr32.seurat <- CreateSeuratObject(counts = bmr32.10x, project = "BMR-32", min.features = 2000, min.cells = 10 )  
bms32.seurat <- CreateSeuratObject(counts = bms32.10x, project = "BMS-32", min.features = 2000, min.cells = 10 )  

# Merge datasets into one single seurat object
bm32.minfeature2000.seurat <- merge(bmy32.seurat, c(bmr32.seurat, bms32.seurat), 
                 add.cell.ids = c("BMY-32", "BMR-32", "BMS-32"))

table(bm32.minfeature2000.seurat$orig.ident)
# min feature = 500
# BMR-32 BMS-32 BMY-32 
# 16183   9095   6981 

# min feature = 2000
# BMR-32 BMS-32 BMY-32 
# 12355   8688   4263 

dim(bm32.minfeature2000.seurat)
# min.feature=500: 21998 32259
# min.feature=2000: 21799 25306

rm(bmy32.10x,bmr32.10x,bms32.10x)
rm(bmy32.seurat,bmr32.seurat,bms32.seurat)
gc()

## Normalization round 1
bm32.minfeature2000.seurat <- SCTransform(bm32.minfeature2000.seurat, assay = 'RNA')

## PCA round 1
bm32.minfeature2000.seurat <- RunPCA(object = bm32.minfeature2000.seurat ) 


mat <- Seurat::GetAssayData(bm32.minfeature2000.seurat, assay = "SCT", slot = "scale.data")
pca <- bm32.minfeature2000.seurat[["pca"]]
total_variance <- sum(matrixStats::rowVars(mat))   # Get the total variance:
eigValues = (pca@stdev)^2  ## EigenValues
varExplained = eigValues / total_variance

plot(c(1:50),varExplained)

## UMAP   round 1
bm32.minfeature2000.seurat <- FindNeighbors(object = bm32.minfeature2000.seurat)
bm32.minfeature2000.seurat <- FindClusters(object = bm32.minfeature2000.seurat, resolution = 0.3)
bm32.minfeature2000.seurat <- RunUMAP(object = bm32.minfeature2000.seurat, dims = 1:20)

##### filter and re-normalize/cluster/etc
# Remove cluster 12  (poor nuclear lysis )
bm32.minfeature2000.seurat <- subset(bm32.minfeature2000.seurat, idents = c('12'), invert=TRUE)

# How many cells are in each sample
table(bm32.minfeature2000.seurat$orig.ident)

# BMR-32 BMS-32 BMY-32 
# 12307   8688   3947 

# How many cells are in each cluster
table(Idents(bm32.minfeature2000.seurat))

## Normalization after filter
bm32.minfeature2000.seurat <- SCTransform(bm32.minfeature2000.seurat, assay = 'RNA')

## PCA after filter
bm32.minfeature2000.seurat <- RunPCA(object = bm32.minfeature2000.seurat ) 


mat <- Seurat::GetAssayData(bm32.minfeature2000.seurat, assay = "SCT", slot = "scale.data")
pca <- bm32.minfeature2000.seurat[["pca"]]
total_variance <- sum(matrixStats::rowVars(mat))   # Get the total variance:
eigValues = (pca@stdev)^2  ## EigenValues
varExplained = eigValues / total_variance

plot(c(1:50),varExplained)

## UMAP
bm32.minfeature2000.seurat <- FindNeighbors(object = bm32.minfeature2000.seurat)
bm32.minfeature2000.seurat <- FindClusters(object = bm32.minfeature2000.seurat, resolution = 0.11)
bm32.minfeature2000.seurat <- RunUMAP(object = bm32.minfeature2000.seurat, dims = 1:20)

saveRDS(bm32.minfeature2000.seurat, file="~/Dropbox (UMass Medical School)/Zinger_at_UMass/Silvia Coreva/Tammy/singleCell/BM32_SingleCell/Iteration3_SeuratMinFeature/bm32.minfeature2000.seurat.RDS")




## Plot PCA, UMAP
DimPlot(bm32.minfeature2000.seurat, reduction = "pca" )
DimPlot(bm32.minfeature2000.seurat, reduction = "pca" , group.by = "orig.ident", cols=c("#d62728","#1f77b4","#e7ba52") )

DimPlot(bm32.minfeature2000.seurat, reduction = "umap", label=T )
DimPlot(bm32.minfeature2000.seurat, reduction = "umap", group.by = "orig.ident", cols=c("#d62728","#1f77b4","#e7ba52") )
 
## Plot gene expression on UMAP projection
FeaturePlot(bm32.minfeature2000.seurat, features = c("CEBPA"), ncol = 1, 
            cols=c("#eeeeee","#d62728"),  
            pt.size = 0.5, order=T) & NoAxes()

## Export UMAP coordinates for Loupe browser
bm32.minfeature2000.seurat.umap.coordinates <- as.data.frame(Embeddings(object = bm32.minfeature2000.seurat[["umap"]])) %>%
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

write.csv(bm32.minfeature2000.seurat.umap.coordinates, file="~/Documents/GitHub/TammyNguyen_SingleCell/BM-32/Iteration 3/bm32.minfeature2000.seurat.umap.csv", row.names=F, quote=F)
  
  

# resolution
# 0.1: 5 clusters #
# 0.1075: 6 clusters #
# 0.11: 7 clusters #
# 0.2: 8 clusters #
# 0.36: 9 clusters # 
# 0.38: 10 clusters #  
# 0.4: 12 clusters # 
# 0.5: 14 clusters #
 

bm32.minfeature2000.seurat.clusters <- bm32.minfeature2000.seurat@meta.data %>% 
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
  select(barcode,`SCT_snn_res.0.1`, `SCT_snn_res.0.1075`,`SCT_snn_res.0.11`,
         `SCT_snn_res.0.2`,`SCT_snn_res.0.36`,`SCT_snn_res.0.38`,`SCT_snn_res.0.4`,`SCT_snn_res.0.5`)


colnames(bm32.minfeature2000.seurat.clusters) <-
  c("barcode","Seurat_min.feature2000_5clusters","Seurat_min.feature2000_6clusters",
    "Seurat_min.feature2000_7clusters","Seurat_min.feature2000_8clusters",
    "Seurat_min.feature2000_9clusters","Seurat_min.feature2000_10clusters",
    "Seurat_min.feature2000_12clusters","Seurat_min.feature2000_14clusters")

write.csv(bm32.minfeature2000.seurat.clusters, file="~/Documents/GitHub/TammyNguyen_SingleCell/BM-32/Iteration 3/bm32.merge.minfeature2000.seurat.clusters.csv", row.names=F, quote=F)

 

##----- QC -----
# https://nbisweden.github.io/workshop-scRNAseq/labs/compiled/seurat/seurat_01_qc.html
bm32.minfeature2000.seurat <- PercentageFeatureSet(bm32.minfeature2000.seurat, "^MT-", col.name = "percent_mito")
bm32.minfeature2000.seurat <- PercentageFeatureSet(bm32.minfeature2000.seurat, "^RP[SL]", col.name = "percent_ribo")
bm32.minfeature2000.seurat <- PercentageFeatureSet(bm32.minfeature2000.seurat, "^HB[^(P)]", col.name = "percent_hemoglobin")
# hemoglobin genes - includes all genes starting with HB except HBP.


# Violin plot by cluster
VlnPlot(bm32.minfeature2000.seurat,  
        features = c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo" ),  
        pt.size = 0, ncol = 4) + 
  NoLegend()

# Violin plot by sample
VlnPlot(bm32.minfeature2000.seurat, group.by = "orig.ident", 
        features = c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo" ),  
        pt.size = 0, ncol = 4) + 
  NoLegend()


##----- Differential Expression ------
# Cluster markers 
# (ident.1 cluster vs all others, genes must be detected in at least 25% of cells in one of the comparison group)
# (log2 fold change >1 or <-1 results are returned)
marker.clust9 <- FindMarkers(bm32.minfeature2000.seurat, ident.1 = "9", min.pct = 0.25, logfc.threshold = 1)
marker.clust8 <- FindMarkers(bm32.minfeature2000.seurat, ident.1 = "8", min.pct = 0.25, logfc.threshold = 1) 
marker.clust7 <- FindMarkers(bm32.minfeature2000.seurat, ident.1 = "7", min.pct = 0.25, logfc.threshold = 1) 

# Comparing two clusters 
marker.clust0v5 <- FindMarkers(bm32.minfeature2000.seurat, ident.1 = "0",ident.2 = "5", min.pct = 0.25, logfc.threshold = 1)
