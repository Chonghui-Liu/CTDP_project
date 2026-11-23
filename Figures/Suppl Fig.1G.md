### The Chi-square algorithm is executed
```r
      load("/data1/liuchonghui/CTDP/test/GSE182227_sc_data.rdata")
      
      chi_square_test <- function(df) {
        
        if (!all(c("seurat_clusters", "source_phen") %in% colnames(df))) {
          stop("data_frame must contain 'cluster' and 'source_phen' column。")
        }
        
        clusters <- unique(df$seurat_clusters)
        
        p_values <- sapply(clusters, function(cl) {
          
          in_cluster <- df$seurat_clusters == cl
          
          contingency_table <- table(
            Cluster = factor(in_cluster, levels = c(FALSE, TRUE)),
            source_phen = factor(df$source_phen, levels = c(0, 1))
          )
          
          test_result <- suppressWarnings(
            chisq.test(contingency_table, correct = FALSE)
          )
          
          return(test_result$p.value)
        })
        
        FDR <- p.adjust(p_values, method = "BH")
        
        positive_rate <- sapply(clusters, function(cl) {
          mean(df$source_phen[df$cluster == cl] == 1, na.rm = TRUE)
        })
        
        result <- data.frame(
          Cluster = clusters,
          Positive_Rate = positive_rate,
          P_Value = as.vector(p_values),
          FDR = as.vector(FDR),
          row.names = paste0("cluster", clusters)
        )
        
        return(result)
      }
      
      
      library(Seurat)
      library(dplyr)
      
      meta_data <- GSE182227_sc_data@meta.data
      
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
      
      GSE182227_chi_square_result <- chi_square_test(df)
      
      save(
        GSE182227_chi_square_result,
        file = "/data1/liuchonghui/CTDP/test/GSE182227_chi_square_result.rdata"
      )
      
``````
### Generate an UMAP plot
```r
      load("/data1/liuchonghui/CTDP/test/GSE182227_sc_data.rdata")
      seurat_obj <- GSE182227_sc_data

      load("/data1/liuchonghui/CTDP/test/GSE182227_chi_square_result.rdata")
      
      chi_square_result <- GSE182227_chi_square_result
      
      significant_clusters <- rownames(
        chi_square_result[chi_square_result$FDR < 0.05, ]
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
      
      colnames(pred_labels) <- "chi_square_pred_labels"
      rownames(pred_labels) <- rownames(seurat_obj@meta.data)
      
      seurat_obj <- AddMetaData(seurat_obj, metadata = pred_labels)
      colnames(seurat_obj@meta.data)

      library(Seurat)
      library(ggpubr)
      library(ggplot2)
      cairo_pdf(
        "/data1/liuchonghui/CTDP/figures of graduation/GSE182227_chi_square_pred_UMAP.pdf"
      )
      DimPlot(
        seurat_obj,
        group.by = "chi_square_pred_labels",
        reduction = "umap",
        label = F,
        repel = TRUE,
        pt.size = 0.5
      ) + NoLegend() +
        theme(aspect.ratio = 1, legend.position = "right")
      dev.off()

``````
