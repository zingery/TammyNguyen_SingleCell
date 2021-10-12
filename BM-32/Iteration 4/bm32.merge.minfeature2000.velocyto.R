library(tidyverse)
library(Seurat)
library(SeuratWrappers)
library(velocyto.R)
library(ggplot2)
library(ggrepel)

#------------------------------------------------------------------------------------------------

# Original Seurat object
setwd("BM32_SingleCell/")
bm32.minfeature2000.seurat <- readRDS("Iteration3_SeuratMinFeature/bm32.minfeature2000.seurat.RDS")

# Velocity read
bmr32.velo <- ReadVelocity("~/Dropbox (UMass Medical School)/Zinger_at_UMass/Silvia Coreva/Tammy/singleCell/BM32_SingleCell/Iteration4_Velocyto/BMR-32.loom")
bmy32.velo <- ReadVelocity("~/Dropbox (UMass Medical School)/Zinger_at_UMass/Silvia Coreva/Tammy/singleCell/BM32_SingleCell/Iteration4_Velocyto/BMY-32.loom")
bms32.velo <- ReadVelocity("~/Dropbox (UMass Medical School)/Zinger_at_UMass/Silvia Coreva/Tammy/singleCell/BM32_SingleCell/Iteration4_Velocyto/BMS-32.loom")

bmr32.velo.Seurat <- as.Seurat(x = bmr32.velo) %>% subset(subset = nFeature_spliced > 1000)
bmy32.velo.Seurat <- as.Seurat(x = bmy32.velo) %>% subset(subset = nFeature_spliced > 1000)
bms32.velo.Seurat <- as.Seurat(x = bms32.velo) %>% subset(subset = nFeature_spliced > 1000)

rm(bmr32.velo,bmy32.velo,bms32.velo)

# the cell barcodes need to be reformatted 
bmr32.velo.Seurat <- RenameCells(
  object = bmr32.velo.Seurat,
  new.names = colnames(bmr32.velo.Seurat) %>% str_replace(":","_") %>% str_replace("x","-1")
)
bmy32.velo.Seurat <- RenameCells(
  object = bmy32.velo.Seurat,
  new.names = colnames(bmy32.velo.Seurat) %>% str_replace(":","_") %>% str_replace("x","-1")
)
bms32.velo.Seurat <- RenameCells(
  object = bms32.velo.Seurat,
  new.names = colnames(bms32.velo.Seurat) %>% str_replace(":","_") %>% str_replace("x","-1")
)

# update the orig.ident
bmr32.velo.Seurat@meta.data$orig.ident = colnames(bmr32.velo.Seurat) %>% str_sub(1,6)
bmy32.velo.Seurat@meta.data$orig.ident = colnames(bmy32.velo.Seurat) %>% str_sub(1,6)
bms32.velo.Seurat@meta.data$orig.ident = colnames(bms32.velo.Seurat) %>% str_sub(1,6)

# remove the cells not in the original minfeature2000 object
bmr32.velo.Seurat <- bmr32.velo.Seurat[,colnames(bmr32.velo.Seurat) %in% colnames(bm32.minfeature2000.seurat)]
bmy32.velo.Seurat <- bmy32.velo.Seurat[,colnames(bmy32.velo.Seurat) %in% colnames(bm32.minfeature2000.seurat)]
bms32.velo.Seurat <- bms32.velo.Seurat[,colnames(bms32.velo.Seurat) %in% colnames(bm32.minfeature2000.seurat)]

# need pca for the velocity
bmr32.velo.Seurat <- SCTransform(bmr32.velo.Seurat,assay="spliced")
bmy32.velo.Seurat <- SCTransform(bmy32.velo.Seurat,assay="spliced")
bms32.velo.Seurat <- SCTransform(bms32.velo.Seurat,assay="spliced")

bmr32.velo.Seurat <- RunPCA(object = bmr32.velo.Seurat )
bmy32.velo.Seurat <- RunPCA(object = bmy32.velo.Seurat )
bms32.velo.Seurat <- RunPCA(object = bms32.velo.Seurat )


mat <- Seurat::GetAssayData(bms32.velo.Seurat, assay = "SCT", slot = "scale.data")
pca <- bmy32.velo.Seurat[["pca"]]
total_variance <- sum(matrixStats::rowVars(mat))   # Get the total variance:
eigValues = (pca@stdev)^2  ## EigenValues
varExplained = eigValues / total_variance
plot(c(1:50),varExplained)
rm(mat,pca)

bmr32.velo.Seurat <- FindNeighbors(object = bmr32.velo.Seurat, dims = 1:15)
bmr32.velo.Seurat <- FindClusters(object = bmr32.velo.Seurat, resolution = 0.2)
bmr32.velo.Seurat <- RunUMAP(object = bmr32.velo.Seurat, dims = 1:15)

bmy32.velo.Seurat <- FindNeighbors(object = bmy32.velo.Seurat, dims = 1:15)
bmy32.velo.Seurat <- FindClusters(object = bmy32.velo.Seurat, resolution = 0.2)
bmy32.velo.Seurat <- RunUMAP(object = bmy32.velo.Seurat, dims = 1:15)

bms32.velo.Seurat <- FindNeighbors(object = bms32.velo.Seurat, dims = 1:15)
bms32.velo.Seurat <- FindClusters(object = bms32.velo.Seurat, resolution = 0.2)
bms32.velo.Seurat <- RunUMAP(object = bms32.velo.Seurat, dims = 1:15)



bmr32.velo.Seurat <- RunVelocity(object = bmr32.velo.Seurat, kCells = 25,fit.quantile = 0.02)
bmy32.velo.Seurat <- RunVelocity(object = bmy32.velo.Seurat, kCells = 25,fit.quantile = 0.02)
bms32.velo.Seurat <- RunVelocity(object = bms32.velo.Seurat, kCells = 25,fit.quantile = 0.02)



d3.scale.category20_2 <- c(
  "#1f77b4",
  "#aec7e8",
  "#ff7f0e",
  "#ffbb78",
  "#2ca02c",
  "#98df8a",
  "#d62728",
  "#ff9896",
  "#9467bd",
  "#c5b0d5",
  "#8c564b",
  "#c49c94",
  "#e377c2",
  "#f7b6d2",
  "#7f7f7f",
  "#c7c7c7",
  "#bcbd22",
  "#dbdb8d",
  "#17becf",
  "#9edae5")

ident.colors <- d3.scale.category20_2

names(x = ident.colors) <- levels(x = bmr32.velo.Seurat)
cell.colors <- ident.colors[Idents(object = bmr32.velo.Seurat)]
names(x = cell.colors) <- colnames(x = bmr32.velo.Seurat)

show.velocity.on.embedding.cor(emb = Embeddings(object = bmr32.velo.Seurat, reduction = "umap"), 
                                      vel = Tool(object = bmr32.velo.Seurat, slot = "RunVelocity"), 
                                      n = 250, scale = "sqrt", 
                                      cell.colors = ac(x = cell.colors, alpha = 0.5), 
                                      cex = 0.8, arrow.scale = 4, show.grid.flow = TRUE, 
                                      min.grid.cell.mass = 0.5, grid.n = 50, arrow.lwd = 1, 
                                      do.par = FALSE, cell.border.alpha = 0.1)


names(x = ident.colors) <- levels(x = bmy32.velo.Seurat)
cell.colors <- ident.colors[Idents(object = bmy32.velo.Seurat)]
names(x = cell.colors) <- colnames(x = bmy32.velo.Seurat)

show.velocity.on.embedding.cor(emb = Embeddings(object = bmy32.velo.Seurat, reduction = "umap"), 
                               vel = Tool(object = bmy32.velo.Seurat, slot = "RunVelocity"), 
                               n = 250, scale = "sqrt", 
                               cell.colors = ac(x = cell.colors, alpha = 0.5), 
                               cex = 0.8, arrow.scale = 4, show.grid.flow = TRUE, 
                               min.grid.cell.mass = 0.5, grid.n = 25, arrow.lwd = 1, 
                               do.par = FALSE, cell.border.alpha = 0.1)


names(x = ident.colors) <- levels(x = bms32.velo.Seurat)
cell.colors <- ident.colors[Idents(object = bms32.velo.Seurat)]
names(x = cell.colors) <- colnames(x = bms32.velo.Seurat)

show.velocity.on.embedding.cor(emb = Embeddings(object = bms32.velo.Seurat, reduction = "umap"), 
                               vel = Tool(object = bms32.velo.Seurat, slot = "RunVelocity"), 
                               n = 250, scale = "sqrt", n.cores=4,
                               cell.colors = ac(x = cell.colors, alpha = 0.5), 
                               cex = 0.8, arrow.scale = 4, show.grid.flow = TRUE, 
                               min.grid.cell.mass = 0.5, grid.n = 50, arrow.lwd = 1, 
                               do.par = FALSE, cell.border.alpha = 0.1)




##----- BMY + BMR -----
# merge the 2 velocity objects
bmyr32.velo.Seurat <- merge(bmr32.velo.Seurat, y=c(bmy32.velo.Seurat)) 
rm(bmr32.velo.Seurat,bmy32.velo.Seurat)

bmyr32.velo.Seurat <- SCTransform(bmyr32.velo.Seurat,assay="spliced") 
bmyr32.velo.Seurat <- RunPCA(object = bmyr32.velo.Seurat ) 


mat <- Seurat::GetAssayData(bmyr32.velo.Seurat, assay = "SCT", slot = "scale.data")
pca <- bmyr32.velo.Seurat[["pca"]]
total_variance <- sum(matrixStats::rowVars(mat))   # Get the total variance:
eigValues = (pca@stdev)^2  ## EigenValues
varExplained = eigValues / total_variance
plot(c(1:50),varExplained)
rm(mat,pca)

bmyr32.velo.Seurat <- FindNeighbors(object = bmyr32.velo.Seurat, dims = 1:15)
bmyr32.velo.Seurat <- FindClusters(object = bmyr32.velo.Seurat, resolution = 0.2)
bmyr32.velo.Seurat <- RunUMAP(object = bmyr32.velo.Seurat, dims = 1:15)
bmyr32.velo.Seurat <- RunVelocity(object = bmyr32.velo.Seurat, kCells = 25,fit.quantile = 0.02)

 

names(x = ident.colors) <- levels(x = bmyr32.velo.Seurat)
cell.colors <- ident.colors[Idents(object = bmyr32.velo.Seurat)]
names(x = cell.colors) <- colnames(x = bmyr32.velo.Seurat)

show.velocity.on.embedding.cor(emb = Embeddings(object = bmyr32.velo.Seurat, reduction = "umap"), 
                               vel = Tool(object = bmyr32.velo.Seurat, slot = "RunVelocity"), 
                               n = 250, scale = "sqrt", 
                               cell.colors = ac(x = cell.colors, alpha = 0.5), 
                               cex = 0.8, arrow.scale = 4, show.grid.flow = TRUE, 
                               min.grid.cell.mass = 0.5, grid.n = 50, arrow.lwd = 1, 
                               do.par = FALSE, cell.border.alpha = 0.1)



# UMAP, color by sample (original identifier) 
DimPlot(bmyr32.velo.Seurat, reduction = "umap", 
        group.by = "orig.ident", cols=c("#d62728","#e7ba52"),   pt.size = 0.5 )

FeaturePlot(bmyr32.velo.Seurat, features = c("LEPR","TAGLN","MKI67","CADM3"), ncol = 2, 
            cols=c("#DEDEDE","#d62728"), 
            pt.size = 0.5, order=T) 
 
 
#