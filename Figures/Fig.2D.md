## Fig2.D
### The DA-seq algorithm is executed
```r
  library(Seurat)
  library(DAseq)
  load("/data/chliu/Supervised_Phen_subpopulations_project/test/PBMC_Travaglini_SC_cluster0to9labeled_list.rdata")
  
  da_cells_list <- lapply(PBMC_Travaglini_SC_cluster0to9labeled_list, function(sc_data) {
    da_cells <- getDAcells(
      X = sc_data@reductions$pca@cell.embeddings,
      cell.labels = as.character(sc_data@meta.data$source_phen),
      labels.1 = "1",
      labels.2 = "0",
      k.vector = seq(50, 500, 50)
    )
    return(da_cells)
  })
  save(da_cells_list,file = "/data/chliu/Supervised_Phen_subpopulations_project/test/da_cells_list.rdata")
``````
### Generate an UMAP plot
```r
  load("/data/chliu/Supervised_Phen_subpopulations_project/test/da_cells_list.rdata")
  da_cells <- da_cells_list[[i]]  
  pred_labels <- rep(0, nrow(seurat_obj@meta.data))
  pred_labels[da_cells$da.down] <- 1
  pred_labels <- as.data.frame(pred_labels)
  table(pred_labels)
  colnames(pred_labels) <- "DAseq_pred_labels"
  rownames(pred_labels) <- rownames(seurat_obj@meta.data)
  print(table(pred_labels))
  seurat_obj <- AddMetaData(seurat_obj, metadata = pred_labels)
  seurat_obj@meta.data[1:5,16:ncol(seurat_obj@meta.data)]

  cairo_pdf(paste0("/data/chliu/Supervised_Phen_subpopulations_project/figures of graduation/",names(PBMC_Travaglini_SC_cluster0to9labeled_list)[i],"_DAseq_pred_UMAP.pdf"))
  DimPlot(seurat_obj, group.by = "DAseq_pred_labels", reduction = "umap", label = F, repel = TRUE,pt.size = 0.5) + NoLegend()+theme(aspect.ratio=1, legend.position = "right")
  dev.off() 
``````