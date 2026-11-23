### The scDist algorithm is executed
```r
      library(Seurat)
      library(scDist)
      
      load("/data1/liuchonghui/CTDP/test/GSE182227_sc_data.rdata")
      
      meta_data <- GSE182227_sc_data@meta.data
      table(meta_data$Sequenced_sample)
      cell.labels <- ifelse(
        meta_data$Sequenced_sample %in% c(
          "OP33",
          "OP34",
          "OP35BOT"
        ),
        "1",
        "0"
      )
      
      table(cell.labels)
      GSE182227_sc_data <- AddMetaData(GSE182227_sc_data, metadata = cell.labels,col.name = "sample_phenotype")

      GSE182227_scDist_result <- scDist(
        LayerData(GSE182227_sc_data, assay = "RNA", layer = "scale.data"),
        GSE182227_sc_data@meta.data,
        fixed.effects = "sample_phenotype",
        random.effects = "Patient",
        clusters = "seurat_clusters"
      )
    
      save(
        GSE182227_scDist_result,
        file = "/data1/liuchonghui/CTDP/test/GSE182227_scDist_result.rdata"
      )
``````
### Generate an UMAP plot
```r
      load("/data1/liuchonghui/CTDP/test/GSE182227_sc_data.rdata")
      seurat_obj <- GSE182227_sc_data

      load("/data1/liuchonghui/CTDP/test/GSE182227_scDist_result.rdata")
      
      scDist_result <- GSE182227_scDist_result
      
      significant_clusters <- rownames(
        scDist_result$results[scDist_result$results$p.val < 0.05, ]
      )
      print(significant_clusters)
      
      pred_labels <- as.data.frame(
        ifelse(
          seurat_obj@meta.data$seurat_clusters %in%
            as.numeric(significant_clusters),
          1, 0
        )
      )
      table(pred_labels)
      
      colnames(pred_labels) <- "scDist_pred_labels"
      rownames(pred_labels) <- rownames(seurat_obj@meta.data)
      
      seurat_obj <- AddMetaData(seurat_obj, metadata = pred_labels)
      colnames(seurat_obj@meta.data)

      library(Seurat)
      library(ggpubr)
      library(ggplot2)
      cairo_pdf(
        "/data1/liuchonghui/CTDP/figures of graduation/GSE182227_scDist_pred_UMAP.pdf"
      )
      DimPlot(
        seurat_obj,
        group.by = "scDist_pred_labels",
        reduction = "umap",
        label = F,
        repel = TRUE,
        pt.size = 0.5
      ) + NoLegend() +
        theme(aspect.ratio = 1, legend.position = "right")
      dev.off()
``````
