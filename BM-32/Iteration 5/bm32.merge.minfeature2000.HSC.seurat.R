library(tidyverse) 
library(data.table)
library(Seurat)

bm32.minfeature2000.seurat <- readRDS("~/Dropbox (UMass Medical School)/Zinger_at_UMass/Silvia Coreva/Tammy/singleCell/BM32_SingleCell/Iteration3_SeuratMinFeature/bm32.minfeature2000.seurat.RDS")

gene_info <- fread("~/gene_info.txt")

## ----- CD34+ cells do not form specific spatial cluster in UMAP projection -----
# THY1 = CD90, PTPRC = CD45
FeaturePlot(bm32.minfeature2000.seurat, features = c("THY1","CD34","PTPRC"),  
            cols=c("#cccccc","#d62728"),  
            pt.size = 0.5, order=T) & NoAxes()

FeaturePlot(bm32.minfeature2000.seurat, features = c("CD34"),  
            cols=c("#cccccc","#d62728"),  
            pt.size = 0.5, order=T) & NoAxes()

# CD34+ cells are visualized as outliers                                  
VlnPlot(bm32.minfeature2000.seurat, features = c("CD34"), group.by = "orig.ident", slot = "scale.data")                                      

## ----- There are 244 cells with any detection of CD34 at raw count -----
bm32.minfeature2000.cd34 <- GetAssayData(object = bm32.minfeature2000.seurat, slot = "counts")['CD34',]
bm32.minfeature2000.cd34.pos.cells <- bm32.minfeature2000.cd34[bm32.minfeature2000.cd34>0]

bm32.minfeature2000.cd34.pos <- subset(bm32.minfeature2000.seurat, cells = names(bm32.minfeature2000.cd34.pos.cells))
VlnPlot(bm32.minfeature2000.cd34.pos, features = c("CD34"), group.by = "orig.ident", slot = "scale.data")                                      

## Breakdown of the CD34+ cells 
bm32.minfeature2000.seurat@meta.data %>% as.data.frame() %>% group_by(orig.ident) %>% summarize(cell_count = n()) %>% ungroup()
#  BMR-32          12307
#  BMS-32           8688
#  BMY-32           3947
bm32.minfeature2000.cd34.pos@meta.data %>% as.data.frame() %>% group_by(orig.ident) %>% summarize(cell_count = n()) %>% ungroup()
#  BMR-32             59  (0.48%)
#  BMS-32            148  (1.70%)
#  BMY-32             37  (0.94%)


bm32.minfeature2000.cd34.pos.norm_count <- GetAssayData(bm32.minfeature2000.cd34.pos, slot = "data")
bm32.minfeature2000.cd34.pos.norm_count <- bm32.minfeature2000.cd34.pos.norm_count %>%
  as.data.frame()

bm32.minfeature2000.cd34.pos.norm_count <- bm32.minfeature2000.cd34.pos.norm_count[rowSums(bm32.minfeature2000.cd34.pos.norm_count[])>1,]

# all-by-all correlation of normalized ocunt
bm32.minfeature2000.cd34.pos.norm_count.cor <- cor(t(bm32.minfeature2000.cd34.pos.norm_count))

# Top correlating genes to CD34 regardless of origin
bm32.minfeature2000.cd34.pos.norm_count.cd34.top12cor <-bm32.minfeature2000.cd34.pos.norm_count.cor["CD34",] %>% sort(decreasing = T) %>% head(12)
FeatureScatter(bm32.minfeature2000.cd34.pos, feature1 = "CD34",feature2="CLDN5")


FeaturePlot(bm32.minfeature2000.seurat, features = names(bm32.minfeature2000.cd34.pos.norm_count.cd34.top10cor),
            cols=c("#cccccc","#d62728"),  
            pt.size = 0.5, order=T) & NoAxes()

data.frame(as.list(bm32.minfeature2000.cd34.pos.norm_count.cd34.top12cor)) %>% t() %>% as.data.frame()%>% tibble::rownames_to_column(var="hgnc_symbol") %>% left_join(gene_info) %>% 
  separate(description, into = c("short_desc","blah"), sep="\\[") %>% select(hgnc_symbol, short_desc) %>% View()


## CD34+ THY1+
bm32.minfeature2000.cd34.pos.thy1 <- GetAssayData(object = bm32.minfeature2000.cd34.pos, slot = "counts")['THY1',]
bm32.minfeature2000.cd34.pos.thy1.cells <- bm32.minfeature2000.cd34.pos.thy1[bm32.minfeature2000.cd34.pos.thy1>0]

bm32.minfeature2000.cd34.pos.thy1.pos <- subset(bm32.minfeature2000.cd34.pos, cells = names(bm32.minfeature2000.cd34.pos.thy1.cells))

bm32.minfeature2000.cd34.pos.thy1.pos@meta.data %>% as.data.frame() %>% group_by(orig.ident) %>% summarize(cell_count = n()) %>% ungroup()

# BMR-32             41
# BMS-32            123
# BMY-32             26
 

## ----- Differential expression analysis between CD34+ cells of different origins -----
de.bmy <- FindMarkers(bm32.minfeature2000.cd34.pos, ident.1 = "BMY-32", group.by = "orig.ident")  
de.bmr <- FindMarkers(bm32.minfeature2000.cd34.pos, ident.1 = "BMR-32", group.by = "orig.ident")  
de.bms <- FindMarkers(bm32.minfeature2000.cd34.pos, ident.1 = "BMS-32", group.by = "orig.ident")  

de.bmy_bmr <- FindMarkers(bm32.minfeature2000.cd34.pos, ident.1 = "BMY-32", ident.2 = "BMR-32", group.by = "orig.ident")  

de.bmy.top12 <- de.bmy %>% slice_max(order_by = avg_log2FC, n=12)
de.bmr.top12 <- de.bmr %>% slice_max(order_by = avg_log2FC, n=12)
de.bmybmr.top12 <- de.bms %>% slice_min(order_by = avg_log2FC, n=12)

de.bmy_bmr.top12 <- de.bmy_bmr %>% slice_min(order_by = avg_log2FC, n=12)


VlnPlot(bm32.minfeature2000.cd34.pos, features = rownames(de.bmy.top12), group.by = "orig.ident", slot = "scale.data")                                  
VlnPlot(bm32.minfeature2000.seurat, features = rownames(de.bmy.top12), group.by = "orig.ident", slot = "scale.data", pt.size = 0)                                      

VlnPlot(bm32.minfeature2000.cd34.pos, features = rownames(de.bmr.top12), group.by = "orig.ident", slot = "scale.data")                                      
VlnPlot(bm32.minfeature2000.seurat, features = rownames(de.bmr.top12), group.by = "orig.ident", slot = "scale.data", pt.size = 0)                                      

VlnPlot(bm32.minfeature2000.cd34.pos, features = rownames(de.bmybmr.top12), group.by = "orig.ident", slot = "scale.data")                                  
VlnPlot(bm32.minfeature2000.seurat, features = rownames(de.bmybmr.top12), group.by = "orig.ident", slot = "scale.data", pt.size = 0)                                      


VlnPlot(bm32.minfeature2000.seurat, features = rownames(de.bmy_bmr.top12), group.by = "orig.ident", slot = "scale.data", pt.size = 0)                                      


## ----- Differential expression analysis between CD34+ cells and CD34- of different origins -----

bm32.minfeature2000.cd34.neg.cells <- bm32.minfeature2000.cd34[bm32.minfeature2000.cd34==0]
 
bm32.minfeature2000.cd34.pos.bmy <- bm32.minfeature2000.cd34.pos.cells[grepl("BMY",names(bm32.minfeature2000.cd34.pos.cells))]
bm32.minfeature2000.cd34.neg.bmy <- bm32.minfeature2000.cd34.neg.cells[grepl("BMY",names(bm32.minfeature2000.cd34.neg.cells))]

tmp <- bm32.minfeature2000.seurat@meta.data %>% as.data.frame() %>%
  tibble::rownames_to_column(var="cell_name") %>%
  mutate(cd34 = ifelse( cell_name %in% names(bm32.minfeature2000.cd34.pos.cells),"pos","neg")) %>%
  mutate(orig_cd34 = case_when(
   grepl("BMY",cell_name) ~ paste0("BMY_",cd34),
   grepl("BMR",cell_name) ~ paste0("BMR_",cd34),
   grepl("BMS",cell_name) ~ paste0("BMS_",cd34), 
  ))
bm32.minfeature2000.seurat@meta.data$orig_cd34 = tmp$orig_cd34

de.bmy.cd34 <- FindMarkers(bm32.minfeature2000.seurat, ident.1 = "BMY_pos", ident.2 ="BMY_neg", group.by = "orig_cd34", logfc.threshold = 0.5 )  
de.bmr.cd34 <- FindMarkers(bm32.minfeature2000.seurat, ident.1 = "BMR_pos", ident.2 ="BMR_neg", group.by = "orig_cd34", logfc.threshold = 0.5 )  

write_tsv(de.bmy.cd34 %>% tibble::rownames_to_column(var="gene"), file = "~/Dropbox (UMass Medical School)/Zinger_at_UMass/Silvia Coreva/Tammy/singleCell/BM32_SingleCell/Iteration5_CD34_positive/de.bmy.cd34.tsv" )
write_tsv(de.bmr.cd34  %>% tibble::rownames_to_column(var="gene"), file = "~/Dropbox (UMass Medical School)/Zinger_at_UMass/Silvia Coreva/Tammy/singleCell/BM32_SingleCell/Iteration5_CD34_positive/de.bmr.cd34.tsv" )

de.bmy.cd34.top12 <- de.bmy.cd34 %>% slice_max(order_by = avg_log2FC, n=12)
de.bmr.cd34.top12 <- de.bmr.cd34 %>% slice_max(order_by = avg_log2FC, n=12)

VlnPlot(bm32.minfeature2000.seurat, features = rownames(de.bmy.cd34.top12), 
        group.by = "orig_cd34", slot = "scale.data", pt.size = 0)                                      

VlnPlot(bm32.minfeature2000.seurat, features = rownames(de.bmr.cd34.top12), 
        group.by = "orig_cd34", slot = "scale.data", pt.size = 0)                                      


de.bmy.cd34.top12 %>% tibble::rownames_to_column(var="hgnc_symbol") %>% left_join(gene_info) %>% 
  separate(description, into = c("short_desc","extra"), sep="\\[") %>% select(hgnc_symbol, short_desc) %>% View()

de.bmr.cd34.top12 %>% tibble::rownames_to_column(var="hgnc_symbol") %>% left_join(gene_info) %>% 
  separate(description, into = c("short_desc","extra"), sep="\\[") %>% select(hgnc_symbol, short_desc) %>% View()
