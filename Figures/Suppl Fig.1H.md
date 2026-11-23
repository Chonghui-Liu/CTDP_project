### The Pearson algorithm is executed
```r
      load("/data1/liuchonghui/CTDP/test/GSE182227_sc_data.rdata")
      
      library(Seurat)
      library(dplyr)
      library(caret)
      
      meta_data <- GSE182227_sc_data@meta.data

      pearson_correlation <- function(df) {
        
        df$cluster <- as.factor(df$seurat_clusters)
        df$phenotype <- as.factor(df$source_phen)

        x_cluster <- model.matrix(~ cluster - 1, data = df)
        x <- as.matrix(x_cluster)

        y <- as.numeric(as.character(df$phenotype))

        pearson_coefficients <- numeric(ncol(x))
        p_values <- numeric(ncol(x))
 
        for (j in 1:ncol(x)) {
          
          test_result <- cor.test(
            x[, j],
            y,
            method = "pearson"
          )
          
          pearson_coefficients[j] <- test_result$estimate
          p_values[j] <- test_result$p.value
        }

        FDR <- p.adjust(p_values, method = "BH")

        result <- data.frame(
          Coefficient = pearson_coefficients,
          P_Value = p_values,
          FDR = FDR
        )
        
        rownames(result) <- colnames(x)
        
        return(result)
      }

      df <- data.frame(
        seurat_clusters = meta_data$seurat_clusters,
        source_phen = meta_data$Sequenced_sample
      )
      rownames(df) <- rownames(meta_data)

      table(df$source_phen)
      df$source_phen <- factor(
        ifelse(df$source_phen %in% c("OP33", "OP34", "OP35BOT"), 1, 0),
        levels = c(0, 1)
      )
      table(df$source_phen)

      GSE182227_pearson_result <- pearson_correlation(df)

      save(
        GSE182227_pearson_result,
        file = "/data1/liuchonghui/CTDP/test/GSE182227_pearson_result.rdata"
      )

      print(GSE182227_pearson_result)
``````
### Generate an UMAP plot
```r
      load("/data1/liuchonghui/CTDP/test/GSE182227_sc_data.rdata")
      seurat_obj <- GSE182227_sc_data
      load("/data1/liuchonghui/CTDP/test/GSE182227_pearson_result.rdata")
      
      pearson_result <- GSE182227_pearson_result
      
      significant_clusters <- rownames(
        pearson_result[pearson_result$FDR < 0.05, ]
      )
      print(significant_clusters)
      
      pred_labels <- as.data.frame(
        ifelse(
          seurat_obj@meta.data$seurat_clusters %in%
            as.numeric(gsub("cluster", "", significant_clusters)),
          1, 0
        )
      )
      table(pred_labels)
      
      colnames(pred_labels) <- "pearson_pred_labels"
      rownames(pred_labels) <- rownames(seurat_obj@meta.data)
      
      seurat_obj <- AddMetaData(seurat_obj, metadata = pred_labels)
      colnames(seurat_obj@meta.data)

      library(Seurat)
      library(ggpubr)
      library(ggplot2)
      cairo_pdf(
        "/data1/liuchonghui/CTDP/figures of graduation/GSE182227_pearson_pred_UMAP.pdf"
      )
      DimPlot(
        seurat_obj,
        group.by = "pearson_pred_labels",
        reduction = "umap",
        label = F,
        repel = TRUE,
        pt.size = 0.5
      ) + NoLegend() +
        theme(aspect.ratio = 1, legend.position = "right")
      dev.off()

``````
