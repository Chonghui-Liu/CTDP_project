## Fig.2G
### The Chi-square algorithm is executed
```r
    load("E:/project/Supervised_Phen_subpopulations_project/test/PBMC_Travaglini_SC_cluster0to9labeled_list.rdata")
    library(dplyr)
    
    chi_square_test <- function(df) {
      
      if (!all(c("cluster", "phenotype") %in% colnames(df))) {
        stop("data_frame must contain 'cluster' and 'phenotype' column。")
      }
      
      clusters <- unique(df$cluster)
      
      p_values <- sapply(clusters, function(cl) {
        
        in_cluster <- df$cluster == cl
        
        contingency_table <- table(
          Cluster = factor(in_cluster, levels = c(FALSE, TRUE)),
          Phenotype = factor(df$phenotype, levels = c(0, 1))
        )
        
        test_result <- suppressWarnings(
          chisq.test(contingency_table, correct = FALSE)
        )
        
        return(test_result$p.value)
      })
      
      FDR <- p.adjust(p_values, method = "BH")
      
      positive_rate <- sapply(clusters, function(cl) {
        mean(df$phenotype[df$cluster == cl] == 1, na.rm = TRUE)
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
    
    chi_square_result_list <- lapply(
      PBMC_Travaglini_SC_cluster0to9labeled_list,
      function(seurat_obj) {
        
        df <- data.frame(
          cluster = seurat_obj@meta.data$seurat_clusters,
          phenotype = seurat_obj@meta.data$source_phen
        )
        
        chi_square_test(df)
      }
    )
    
    save(
      chi_square_result_list,
      file = "E:/project/Supervised_Phen_subpopulations_project/test/chi_square_result_list.rdata"
    )
    
``````
### Generate an UMAP plot
```r
    load("E:/project/Supervised_Phen_subpopulations_project/test/PBMC_Travaglini_SC_cluster0to9labeled_list.rdata")
    i = 2
    print(names(PBMC_Travaglini_SC_cluster0to9labeled_list)[i])
    seurat_obj <- PBMC_Travaglini_SC_cluster0to9labeled_list[[i]]
    
    
    load("E:/project/Supervised_Phen_subpopulations_project/test/chi_square_result_list.rdata")
    chi_square_result <- chi_square_result_list[[i]]
    significant_clusters <- rownames(chi_square_result[chi_square_result$FDR < 0.05, ])
    print(significant_clusters)
    
    pred_labels <- as.data.frame(
      ifelse(seurat_obj@meta.data$seurat_clusters %in% as.numeric(gsub("cluster", "", significant_clusters)), 1, 0)
    )
    table(pred_labels)
    colnames(pred_labels) <- "chi_square_pred_labels"
    rownames(pred_labels) <- rownames(seurat_obj@meta.data)
    seurat_obj <- AddMetaData(seurat_obj, metadata = pred_labels)
    seurat_obj@meta.data[1:5, 16:ncol(seurat_obj@meta.data)]
    
    cairo_pdf(paste0("E:/project/Supervised_Phen_subpopulations_project/figures of graduation/", names(PBMC_Travaglini_SC_cluster0to9labeled_list)[i], "_chi_square_pred_UMAP.pdf"))
    DimPlot(seurat_obj, group.by = "chi_square_pred_labels", reduction = "umap", label = F, repel = TRUE, pt.size = 0.5) + NoLegend() + theme(aspect.ratio = 1, legend.position = "right")
    dev.off()
       
``````
