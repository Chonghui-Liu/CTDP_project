## Fig5.A,B
### Download and preprocess the single-cell data of COVID-19
```r
  library(Seurat)
  covid_nbt_main <- readRDS("/data/chliu/Supervised_Phen_subpopulations_project/data/COVID-19/covid_nbt_main.rds")
  dim(covid_nbt_main@meta.data)
  head(covid_nbt_main@meta.data)
  table(covid_nbt_main@meta.data$severity)
  head(Idents(covid_nbt_main))
  
  
  expr_data <- covid_nbt_main@assays$RNA@counts
  dim(expr_data)
  meta_data <- covid_nbt_main@meta.data
  dim(meta_data)
  COVID19_sc_data <- CreateSeuratObject(counts = expr_data, meta.data = meta_data, min.cells = 3, min.features = 200)
  gc()
  print(COVID19_sc_data)
  COVID19_sc_data <- subset(COVID19_sc_data, subset = severity != "control")
  gc()
  print(COVID19_sc_data)
  COVID19_sc_data[["percent.mt"]] <- PercentageFeatureSet(COVID19_sc_data, pattern = "^MT-")
  COVID19_sc_data <- subset(COVID19_sc_data, subset = percent.mt < 20)
  dim(COVID19_sc_data@assays[["RNA"]]@counts)
  COVID19_sc_data <- NormalizeData(COVID19_sc_data, normalization.method = "LogNormalize", scale.factor = 10000)
  COVID19_sc_data <- ScaleData(COVID19_sc_data)
  COVID19_sc_data <- FindVariableFeatures(COVID19_sc_data, selection.method = "vst", nfeatures = 2000)
  COVID19_sc_data <- RunPCA(COVID19_sc_data, features = VariableFeatures(object = COVID19_sc_data))
  COVID19_sc_data <- FindNeighbors(COVID19_sc_data, dims = 1:10)
  COVID19_sc_data <- FindClusters(COVID19_sc_data)
  COVID19_sc_data <- RunUMAP(COVID19_sc_data, dims = 1:10)
  save(COVID19_sc_data,file="/data/chliu/Supervised_Phen_subpopulations_project/test/COVID19_sc_data.rdata")
``````
### Downsampling of single-cell data in COVID-19 cases
```r
  load("/data/chliu/Supervised_Phen_subpopulations_project/test/COVID19_sc_data.rdata")
  table(COVID19_sc_data@meta.data$severity)
  table(COVID19_sc_data@meta.data$seurat_clusters)
  table(COVID19_sc_data@meta.data$celltype)
  library(dplyr)
  library(Seurat)
  meta_data <- COVID19_sc_data@meta.data
  target_count <- sum(meta_data$severity == "critical")
  moderate_counts <- meta_data %>%
    filter(severity == "moderate") %>%
    group_by(celltype) %>%
    summarise(cluster_count = n())
  moderate_counts <- moderate_counts %>%
    mutate(target_sample_size = round(cluster_count / sum(cluster_count) * target_count))
  set.seed(123) 
  sampled_moderate <- data.frame()
  for (i in 1:nrow(moderate_counts)) {
    cluster_id <- moderate_counts$celltype[i]
    sample_size <- moderate_counts$target_sample_size[i]
    
    cells_in_cluster <- meta_data %>%
      filter(severity == "moderate", celltype == cluster_id)
    
    sampled_cells <- if (nrow(cells_in_cluster) <= sample_size) {
      cells_in_cluster
    } else {
      cells_in_cluster[sample(1:nrow(cells_in_cluster), sample_size), ]
    }
    
    sampled_moderate <- rbind(sampled_moderate, sampled_cells)
  }
  dim(sampled_moderate)
  sampled_moderate[1:5,11:17]
  critical <- meta_data %>%
    filter(severity == "critical")
  dim(critical)
  downsampled_meta_data <- bind_rows(sampled_moderate, critical)
  COVID19_sc_data_downsampled <- subset(COVID19_sc_data, cells = rownames(downsampled_meta_data))
  table(downsampled_meta_data$severity)
  dim(COVID19_sc_data_downsampled)
  save(COVID19_sc_data_downsampled,file = "/data/chliu/Supervised_Phen_subpopulations_project/test/COVID19_sc_data_downsampled.rdata")
``````
### UAMP plots of single-cell data for cases of COVID-19
```r
  load("/data/chliu/Supervised_Phen_subpopulations_project/test/COVID19_sc_data_downsampled.rdata")
  library(Seurat)
  library(ggpubr)
  library(ggplot2)
  cairo_pdf(paste0("/data/chliu/Supervised_Phen_subpopulations_project/figures of graduation/COVID19_sc_data_downsampled_UMAP.pdf"))
  DimPlot(COVID19_sc_data_downsampled, group.by ="celltype", reduction = "umap", label = T, repel = TRUE,pt.size = 0.5) + NoLegend()+theme(aspect.ratio=1, legend.position = "right")
  dev.off()
  
  cairo_pdf(paste0("/data/chliu/Supervised_Phen_subpopulations_project/figures of graduation/COVID19_sc_data_downsampled_source_phen_UMAP.pdf"))
  DimPlot(COVID19_sc_data_downsampled, group.by = "severity", reduction = "umap", label = F, repel = TRUE,pt.size = 0.5) + NoLegend()+theme(aspect.ratio=1, legend.position = "right")
  dev.off() 
``````