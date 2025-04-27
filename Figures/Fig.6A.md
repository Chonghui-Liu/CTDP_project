## Fig6.A,B
### Download and preprocess the single-cell data of Liver cirrhosis
```r
  library(Seurat)
  load("tissue.rdata")
  expr_matrix <- tissue@raw.data
  meta_data <- tissue@meta.data
  features <- rownames(expr_matrix)
  head(features)
  features <- gsub("_", "-", features) 
  features <- gsub("\\|", "-", features) 
  rownames(expr_matrix) <- features
  Livercirrhosis_sc_data <- CreateSeuratObject(counts = expr_matrix, meta.data = meta_data, min.cells = 3, min.features = 200)
  gc()
  print(Livercirrhosis_sc_data)
  table(Livercirrhosis_sc_data@meta.data$annotation_lineage)
  table(Livercirrhosis_sc_data@meta.data$condition)
  Livercirrhosis_sc_data[["percent.mt"]] <- PercentageFeatureSet(Livercirrhosis_sc_data, pattern = "^MT-")
  Livercirrhosis_sc_data <- subset(Livercirrhosis_sc_data, subset = percent.mt < 20)
  dim(Livercirrhosis_sc_data@assays[["RNA"]]@counts)
  Livercirrhosis_sc_data <- NormalizeData(Livercirrhosis_sc_data, normalization.method = "LogNormalize", scale.factor = 10000)
  Livercirrhosis_sc_data <- ScaleData(Livercirrhosis_sc_data)
  Livercirrhosis_sc_data <- FindVariableFeatures(Livercirrhosis_sc_data, selection.method = "vst", nfeatures = 2000)
  Livercirrhosis_sc_data <- RunPCA(Livercirrhosis_sc_data, features = VariableFeatures(object = Livercirrhosis_sc_data))
  Livercirrhosis_sc_data <- FindNeighbors(Livercirrhosis_sc_data, dims = 1:10)
  Livercirrhosis_sc_data <- FindClusters(Livercirrhosis_sc_data)
  Livercirrhosis_sc_data <- RunUMAP(Livercirrhosis_sc_data, dims = 1:10)
  save(Livercirrhosis_sc_data,file="Livercirrhosis_sc_data.rdata")
``````
### Downsampling of single-cell data in liver cirrhosis cases
```r
  load("Livercirrhosis_sc_data.rdata")
  table(Livercirrhosis_sc_data@meta.data$condition)
  table(Livercirrhosis_sc_data@meta.data$seurat_clusters)
  table(Livercirrhosis_sc_data@meta.data$annotation_lineage)
  library(dplyr)
  library(Seurat)
  meta_data <- Livercirrhosis_sc_data@meta.data
  target_count <- sum(meta_data$condition == "Cirrhotic")
  Uninjured_counts <- meta_data %>%
    filter(condition == "Uninjured") %>%
    group_by(annotation_lineage) %>%
    summarise(cluster_count = n())
  Uninjured_counts <- Uninjured_counts %>%
    mutate(target_sample_size = round(cluster_count / sum(cluster_count) * target_count))
  set.seed(123) 
  sampled_Uninjured <- data.frame()
  for (i in 1:nrow(Uninjured_counts)) {
    cluster_id <- Uninjured_counts$annotation_lineage[i]
    sample_size <- Uninjured_counts$target_sample_size[i]
    
    cells_in_cluster <- meta_data %>%
      filter(condition == "Uninjured", annotation_lineage == cluster_id)
    
    sampled_cells <- if (nrow(cells_in_cluster) <= sample_size) {
      cells_in_cluster
    } else {
      cells_in_cluster[sample(1:nrow(cells_in_cluster), sample_size), ]
    }
    sampled_Uninjured <- rbind(sampled_Uninjured, sampled_cells)
  }
  dim(sampled_Uninjured)
  sampled_Uninjured[1:5,8:16]
  Cirrhotic <- meta_data %>%
    filter(condition == "Cirrhotic")
  dim(Cirrhotic)
  Cirrhotic[1:5,8:16]
  downsampled_meta_data <- bind_rows(sampled_Uninjured, Cirrhotic)
  dim(downsampled_meta_data)
  Livercirrhosis_sc_data_downsampled <- subset(Livercirrhosis_sc_data, cells = rownames(downsampled_meta_data))
  table(downsampled_meta_data$condition)
  dim(Livercirrhosis_sc_data_downsampled)
  save(Livercirrhosis_sc_data_downsampled,file = "Livercirrhosis_sc_data_downsampled.rdata")
``````
### UAMP plots of single-cell data for cases of liver cirrhosis
```r
  load("Livercirrhosis_sc_data_downsampled.rdata")
  library(Seurat)
  library(ggpubr)
  library(ggplot2)
  cairo_pdf(paste0("figures of graduation/Livercirrhosis_sc_data_downsampled_UMAP.pdf"))
  DimPlot(Livercirrhosis_sc_data_downsampled, group.by ="annotation_lineage", reduction = "umap", label = T, repel = TRUE,pt.size = 0.5) + NoLegend()+theme(aspect.ratio=1, legend.position = "right")
  dev.off()
  
  cairo_pdf(paste0("figures of graduation/Livercirrhosis_sc_data_downsampled_source_phen_UMAP.pdf"))
  DimPlot(Livercirrhosis_sc_data_downsampled, group.by = "condition", reduction = "umap", label = F, repel = TRUE,pt.size = 0.5) + NoLegend()+theme(aspect.ratio=1, legend.position = "right")
  dev.off()
``````
