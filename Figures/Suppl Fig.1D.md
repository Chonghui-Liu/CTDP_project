### The DA-seq algorithm is executed
```r
      library(Seurat)
      library(DAseq)
      
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

      GSE182227_DAseq_result <- getDAcells(
        X = GSE182227_sc_data@reductions$pca@cell.embeddings,
        cell.labels = cell.labels,
        labels.1 = "1",
        labels.2 = "0",
        k.vector = seq(50, 500, 50)
      )

      save(
        GSE182227_DAseq_result,
        file = "/data1/liuchonghui/CTDP/test/GSE182227_DAseq_result.rdata"
      )
``````
### Generate an UMAP plot
```r
      load("/data1/liuchonghui/CTDP/test/GSE182227_sc_data.rdata")
      seurat_obj <- GSE182227_sc_data

      load("/data1/liuchonghui/CTDP/test/GSE182227_DAseq_result.rdata")
      
      da_cells <- GSE182227_DAseq_result
      
      pred_labels <- rep(0, nrow(seurat_obj@meta.data))
      pred_labels[da_cells$da.down] <- 1
      
      pred_labels <- as.data.frame(pred_labels)
      table(pred_labels)
      
      colnames(pred_labels) <- "DAseq_pred_labels"
      rownames(pred_labels) <- rownames(seurat_obj@meta.data)
      
      print(table(pred_labels))
      
      seurat_obj <- AddMetaData(seurat_obj, metadata = pred_labels)
      colnames(seurat_obj@meta.data)
      library(Seurat)
      library(ggpubr)
      library(ggplot2)
      cairo_pdf(
        "/data1/liuchonghui/CTDP/figures of graduation/GSE182227_DAseq_pred_UMAP.pdf"
      )
      DimPlot(
        seurat_obj,
        group.by = "DAseq_pred_labels",
        reduction = "umap",
        label = F,
        repel = TRUE,
        pt.size = 0.5
      ) + NoLegend() +
        theme(aspect.ratio = 1, legend.position = "right")
      dev.off()
``````
