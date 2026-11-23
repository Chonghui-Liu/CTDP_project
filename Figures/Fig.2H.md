## Fig.2H
### The Pearson algorithm is executed
```r
    load("E:/project/Supervised_Phen_subpopulations_project/test/PBMC_Travaglini_SC_cluster0to9labeled_list.rdata")
    
    pearson_correlation <- function(df) {
      
      df <- as.data.frame(df@meta.data)
      
      df$cluster <- as.factor(df$seurat_clusters)
      df$phenotype <- df$source_phen
      
      x_cluster <- model.matrix(~ cluster - 1, data = df)
      x <- as.matrix(x_cluster)
      
      y <- df$phenotype
      
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
    
    pearson_result_list <- list()
    
    for (i in seq_along(
      PBMC_Travaglini_SC_cluster0to9labeled_list
    )) {
      
      seurat_obj <-
        PBMC_Travaglini_SC_cluster0to9labeled_list[[i]]
      
      result <- pearson_correlation(seurat_obj)
      
      name <- names(
        PBMC_Travaglini_SC_cluster0to9labeled_list
      )[i]
      
      pearson_result_list[[name]] <- result
    }
    
    save(
      pearson_result_list,
      file = "E:/project/Supervised_Phen_subpopulations_project/test/pearson_result_list.rdata"
    )
    
``````
### Generate an UMAP plot
```r
    load("E:/project/Supervised_Phen_subpopulations_project/test/PBMC_Travaglini_SC_cluster0to9labeled_list.rdata")
    i = 2
    print(names(PBMC_Travaglini_SC_cluster0to9labeled_list)[i])
    seurat_obj <- PBMC_Travaglini_SC_cluster0to9labeled_list[[i]]
    
    
    load("E:/project/Supervised_Phen_subpopulations_project/test/pearson_result_list.rdata")
    pearson_result <- pearson_result_list[[i]]
    significant_clusters <- rownames(pearson_result[pearson_result$FDR < 0.05, ])
    print(significant_clusters)
    
    pred_labels <- as.data.frame(
      ifelse(seurat_obj@meta.data$seurat_clusters %in% as.numeric(gsub("cluster", "", significant_clusters)), 1, 0)
    )
    table(pred_labels)
    colnames(pred_labels) <- "pearson_pred_labels"
    rownames(pred_labels) <- rownames(seurat_obj@meta.data)
    seurat_obj <- AddMetaData(seurat_obj, metadata = pred_labels)
    seurat_obj@meta.data[1:5, 16:ncol(seurat_obj@meta.data)]
    
    cairo_pdf(paste0("E:/project/Supervised_Phen_subpopulations_project/figures of graduation/", names(PBMC_Travaglini_SC_cluster0to9labeled_list)[i], "_pearson_pred_UMAP.pdf"))
    DimPlot(seurat_obj, group.by = "pearson_pred_labels", reduction = "umap", label = F, repel = TRUE, pt.size = 0.5) + NoLegend() + theme(aspect.ratio = 1, legend.position = "right")
    dev.off()
``````
