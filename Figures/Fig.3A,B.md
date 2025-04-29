## Fig3.A,B
### The single cell data were labeled with simulation
```r
library(Seurat)
  assign_labels <- function(seurat_obj, clusters, ratio) {
    meta_data <- seurat_obj@meta.data
    target_cells <- which(meta_data$seurat_clusters %in% clusters)
    
    # 初始化标签列
    meta_data$source_phen <- NA
    meta_data$groundtruth <- 0
    
    # 为目标簇分配标签
    meta_data$groundtruth[target_cells] <- 1
    meta_data$source_phen[target_cells] <- ifelse(runif(length(target_cells)) < (1 - ratio), 1, 0)
    
    # 为非目标簇分配标签
    non_target_cells <- setdiff(1:nrow(meta_data), target_cells)
    meta_data$source_phen[non_target_cells] <- sample(c(0, 1), length(non_target_cells), replace = TRUE)
    
    seurat_obj@meta.data <- meta_data
    return(seurat_obj)
  }
  load("/data/chliu/Supervised_Phen_subpopulations_project/test/PBMC_Travaglini_SC.rdata")
  
  set.seed(123)
  all_clusters <- unique(PBMC_Travaglini_SC@meta.data$seurat_clusters)
  PBMC_Travaglini_SC_sample3clusters_list <- list()
  for (i in 1:10) {
    random_clusters <- sample(all_clusters, size=3)
    cluster_name <- paste("clusters", paste(random_clusters, collapse = "_"), sep = "_")
    print(cluster_name)
    labeled_obj <- assign_labels(PBMC_Travaglini_SC, clusters = random_clusters, ratio = 0.1)
    PBMC_Travaglini_SC_sample3clusters_list[i] <- labeled_obj
    names(PBMC_Travaglini_SC_sample3clusters_list)[i] <- cluster_name
  }
  save(PBMC_Travaglini_SC_sample3clusters_list, file = "/data/chliu/Supervised_Phen_subpopulations_project/test/PBMC_Travaglini_SC_sample3clusters_list.rdata")
    
  
  
``````
### Generate an UMAP plot
```r
  load("/data/chliu/Supervised_Phen_subpopulations_project/test/PBMC_Travaglini_SC_sample3clusters_list.rdata")
  i = 5
  print(names(PBMC_Travaglini_SC_sample3clusters_list)[i])
  seurat_obj <- PBMC_Travaglini_SC_sample3clusters_list[[i]]
  
  library(Seurat)
  library(ggpubr)
  library(ggplot2)
  cairo_pdf(paste0("/data/chliu/Supervised_Phen_subpopulations_project/figures of graduation/",names(PBMC_Travaglini_SC_sample3clusters_list)[i],"_UMAP.pdf"))
  DimPlot(seurat_obj,  reduction = "umap", label = T, repel = TRUE,pt.size = 0.5) + NoLegend()+theme(aspect.ratio=1, legend.position = "right")
  dev.off()
  
  cairo_pdf(paste0("/data/chliu/Supervised_Phen_subpopulations_project/figures of graduation/",names(PBMC_Travaglini_SC_sample3clusters_list)[i],"_source_phen_UMAP.pdf"))
  DimPlot(seurat_obj, group.by = "source_phen", reduction = "umap", label = F, repel = TRUE,pt.size = 0.5) + NoLegend()+theme(aspect.ratio=1, legend.position = "right")
  dev.off()  
  
  
``````