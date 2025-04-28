## Fig2.F
### The scDist algorithm is executed
```r
  library(Seurat)
  library(scDist)
  split_source_phen <- function(sc_data) {
    if (!"source_phen" %in% colnames(sc_data@meta.data)) {
      stop("source_phen column does not exit in sc_data@meta.data")
    }
    
    source_phen_1_indices <- which(sc_data@meta.data$source_phen == 1)
    source_phen_0_indices <- which(sc_data@meta.data$source_phen == 0)
    
    set.seed(123)  
    group_1_indices <- sample(source_phen_1_indices, length(source_phen_1_indices) / 2)
    group_2_indices <- setdiff(source_phen_1_indices, group_1_indices)
    
    group_3_indices <- sample(source_phen_0_indices, length(source_phen_0_indices) / 2)
    group_4_indices <- setdiff(source_phen_0_indices, group_3_indices)
    
    sc_data@meta.data$patients <- NA  
    sc_data@meta.data$patients[group_1_indices] <- "T1"
    sc_data@meta.data$patients[group_2_indices] <- "T2"
    sc_data@meta.data$patients[group_3_indices] <- "C1"
    sc_data@meta.data$patients[group_4_indices] <- "C2"
    
    return(sc_data)
  }
  
  load("/data/chliu/Supervised_Phen_subpopulations_project/test/PBMC_Travaglini_SC_cluster0to9labeled_list.rdata")
  scDist_result_list <- lapply(PBMC_Travaglini_SC_cluster0to9labeled_list, function(sc_data) {
    sc_data <- split_source_phen(sc_data)
    scDist_result <- scDist(as.matrix(sc_data@assays$RNA@scale.data),
                            sc_data@meta.data,
                            fixed.effects = c("source_phen", "patients"),
                            clusters = "seurat_clusters",
                            d = 8)
    return(scDist_result)
  })
  save(scDist_result_list,file = "/data/chliu/Supervised_Phen_subpopulations_project/test/scDist_result_list.rdata")
   
 
 
``````
### Generate an UMAP plot
```r
  load("/data/chliu/Supervised_Phen_subpopulations_project/test/scDist_result_list.rdata")
  scDist_result <- scDist_result_list[[i]]
  significant_clusters <- rownames(scDist_result$results[scDist_result$results$p.val < 0.05, ])
  print(significant_clusters)
  pred_labels <- as.data.frame(
    ifelse(seurat_obj@meta.data$seurat_clusters %in% as.numeric(significant_clusters), 1, 0)
  )
  table(pred_labels)
  colnames(pred_labels) <- "scDist_pred_labels"
  rownames(pred_labels) <- rownames(seurat_obj@meta.data)
  seurat_obj <- AddMetaData(seurat_obj, metadata = pred_labels)
  seurat_obj@meta.data[1:5,16:ncol(seurat_obj@meta.data)]

  cairo_pdf(paste0("/data/chliu/Supervised_Phen_subpopulations_project/figures of graduation/",names(PBMC_Travaglini_SC_cluster0to9labeled_list)[i],"_scDist_pred_UMAP.pdf"))
  DimPlot(seurat_obj, group.by = "scDist_pred_labels", reduction = "umap", label = F, repel = TRUE,pt.size = 0.5) + NoLegend()+theme(aspect.ratio=1, legend.position = "right")
  dev.off() 
``````