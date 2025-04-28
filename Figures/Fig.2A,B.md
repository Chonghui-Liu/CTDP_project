## Fig2.A,B
### Preprocessing of PBMC_Travaglini single-cell data
```r
  library(Seurat)
  load("D:/GWAS_SC_project/data/SC/PBMC_Travaglini/droplet_normal_lung_blood_seurat_ntiss10x.P1.anno.20191002.RC4.Robj")
  ntiss10x.P1.anno <- UpdateSeuratObject(ntiss10x.P1.anno)
  ntiss10x.P1.anno_blood <- subset(ntiss10x.P1.anno, subset= tissue=="blood")
  gc()
  dim(ntiss10x.P1.anno@assays$RNA@counts)
  dim(ntiss10x.P1.anno_blood@assays$RNA@counts)
  table(ntiss10x.P1.anno_blood$tissue)
  load("D:/GWAS_SC_project/data/SC/PBMC_Travaglini/droplet_normal_lung_blood_seurat_ntiss10x.P3.anno.20191002.RC4.Robj")
  ntiss10x.P3.anno <- UpdateSeuratObject(ntiss10x.P3.anno)
  ntiss10x.P3.anno_blood <- subset(ntiss10x.P3.anno, subset= tissue=="blood")
  gc()
  dim(ntiss10x.P3.anno@assays$RNA@counts)
  dim(ntiss10x.P3.anno_blood@assays$RNA@counts)
  table(ntiss10x.P3.anno_blood$tissue)
  PBMC_Travaglini_Seurat <- merge(ntiss10x.P1.anno_blood,ntiss10x.P3.anno_blood)
  dim(PBMC_Travaglini_Seurat@assays$RNA@counts)
  PBMC_Travaglini_Seurat@assays$RNA@counts[1:5,1:3]
  counts <- PBMC_Travaglini_Seurat@assays$RNA@counts
  dim(counts)
  counts[1:5,1:5]
  PBMC_SC <- CreateSeuratObject(
    counts = counts, 
    meta.data=as.data.frame(PBMC_Travaglini_Seurat@meta.data),
    min.cells = 3, 
    min.features = 200)
  PBMC_SC@meta.data[1:5,]
  dim(PBMC_SC@meta.data)
  table(PBMC_SC@meta.data$free_annotation)
  Idents(PBMC_SC)<-PBMC_SC@meta.data$free_annotation
  table(Idents(PBMC_SC))
  cell_type_counts <- table(Idents(PBMC_SC))
  selected_cell_types <- names(cell_type_counts[cell_type_counts > 10])
  PBMC_SC <- subset(PBMC_SC, idents = selected_cell_types)
  dim(PBMC_SC)
  table(PBMC_SC$free_annotation)
  dim(counts);dim(PBMC_SC@assays[["RNA"]]@counts)
  PBMC_SC[["percent.mt"]] <- PercentageFeatureSet(PBMC_SC, pattern = "^MT-")
  PBMC_SC <- subset(PBMC_SC, subset = percent.mt < 20)
  dim(PBMC_SC@assays[["RNA"]]@counts)
  PBMC_SC <- NormalizeData(PBMC_SC, normalization.method = "LogNormalize", scale.factor = 10000)#这一步运行的结果保存在PBMC_SC@assays$RNA@data
  PBMC_Travaglini_SC <- ScaleData(PBMC_SC)#scPagwas代码把ScaleData函数写前面，NormalizeData函数写后面，这是错误的。#这一步的结果保存在PBMC_SC@assays$RNA@scale.data
  saveRDS(PBMC_Travaglini_SC,file="D:/GWAS_SC_project/test/PBMC_Travaglini_SC.rds")
  
  library(Seurat)
  PBMC_Travaglini_SC <- readRDS("D:/GWAS_SC_project/test/PBMC_Travaglini_SC.rds")
  PBMC_Travaglini_SC <- FindVariableFeatures(PBMC_Travaglini_SC, selection.method = "vst", nfeatures = 2000)
  PBMC_Travaglini_SC <- RunPCA(PBMC_Travaglini_SC, features = VariableFeatures(object = PBMC_Travaglini_SC))
  PBMC_Travaglini_SC <- FindNeighbors(PBMC_Travaglini_SC, dims = 1:10)
  PBMC_Travaglini_SC <- FindClusters(PBMC_Travaglini_SC)
  PBMC_Travaglini_SC <- RunUMAP(PBMC_Travaglini_SC, dims = 1:10)
  save(PBMC_Travaglini_SC,file = "D:/Supervised_Phen_subpopulations_project/test/PBMC_Travaglini_SC.rdata")
  library(ggpubr)
  library(ggplot2)
  pdf("D:/Supervised_Phen_subpopulations_project/figures of graduation/PBMC_Travaglini_UMAP.pdf")
  DimPlot(PBMC_Travaglini_SC, reduction = "umap", label = TRUE, repel = TRUE,pt.size = 0.5) + NoLegend()+theme(aspect.ratio=1)
  dev.off()

``````
### The single cell data were labeled with simulation
```r
  load("D:/Supervised_Phen_subpopulations_project/test/PBMC_Travaglini_SC.rdata")
  assign_labels <- function(seurat_obj, clusters, ratio) {
    set.seed(123)
    meta_data <- seurat_obj@meta.data
    target_cells <- which(meta_data$seurat_clusters %in% clusters)
    
    meta_data$source_phen <- NA
    meta_data$groundtruth <- 0
    
    meta_data$groundtruth[target_cells] <- 1
    meta_data$source_phen[target_cells] <- ifelse(runif(length(target_cells)) < (1 - ratio), 1, 0)
    
    non_target_cells <- setdiff(1:nrow(meta_data), target_cells)
    meta_data$source_phen[non_target_cells] <- sample(c(0, 1), length(non_target_cells), replace = TRUE)
    
    seurat_obj@meta.data <- meta_data
    return(seurat_obj)
  }
  
  PBMC_Travaglini_SC_cluster1labeled <- assign_labels(PBMC_Travaglini_SC, clusters = c(1), ratio = 0.1)
  meta_pos <- PBMC_Travaglini_SC_cluster1labeled@meta.data[WhichCells(PBMC_Travaglini_SC_cluster1labeled, idents = 1), ]
  table(meta_pos$source_phen)
  meta_neg <- PBMC_Travaglini_SC_cluster1labeled@meta.data[WhichCells(PBMC_Travaglini_SC_cluster1labeled, idents = 2), ]
  table(meta_neg$source_phen)
  meta_neg <- PBMC_Travaglini_SC_cluster1labeled@meta.data[WhichCells(PBMC_Travaglini_SC_cluster1labeled, idents = 9), ]
  table(meta_neg$source_phen)
  save(PBMC_Travaglini_SC_cluster1labeled,file = "D:/Supervised_Phen_subpopulations_project/test/PBMC_Travaglini_SC_cluster1labeled.rdata")
    
  
``````
### Generate an UMAP plot
```r
  load("D:/Supervised_Phen_subpopulations_project/test/PBMC_Travaglini_SC_cluster1labeled.rdata")
  library(Seurat)
  library(ggpubr)
  library(ggplot2)
  pdf("D:/Supervised_Phen_subpopulations_project/figures of graduation/PBMC_Travaglini_sourceSamplelabel_UMAP.pdf")
  DimPlot(PBMC_Travaglini_SC_cluster1labeled, group.by = "source_phen", reduction = "umap", label = F, repel = TRUE,pt.size = 0.5) + NoLegend()+theme(aspect.ratio=1, legend.position = "right")
  dev.off()  
  
``````