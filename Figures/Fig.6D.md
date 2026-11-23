## Fig6.D
```r
  load("Livercirrhosis_logistic_result.rdata")
  load("Livercirrhosis_sc_data_downsampled.rdata")
  logistic_result <- Livercirrhosis_logistic_result
  significant_clusters <- rownames(logistic_result[logistic_result$FDR < 0.05, ])
  print(significant_clusters)
  pred_labels <- as.data.frame(
    ifelse(Livercirrhosis_sc_data_downsampled@meta.data$annotation_lineage %in% as.character(gsub("cluster", "", significant_clusters)), 1, 0)
  )
  table(pred_labels)
  colnames(pred_labels) <- "logistic_pred_labels"
  rownames(pred_labels) <- rownames(Livercirrhosis_sc_data_downsampled@meta.data)
  Livercirrhosis_sc_data_downsampled <- AddMetaData(Livercirrhosis_sc_data_downsampled, metadata = pred_labels)
  Livercirrhosis_sc_data_downsampled@meta.data[1:5,8:17]
  
  
  options(future.globals.maxSize = 5 * 1024^3)
  library(future)
  plan("multicore", workers = 50)
  Cirrhotic_associatedcells_sig <- Seurat::FindMarkers(Livercirrhosis_sc_data_downsampled, ident.1 = 1, group.by = "logistic_pred_labels", logfc.threshold=0, min.pct=0,random.seed=123)
  Cirrhotic_associatedcells_sig$names <- rownames(Cirrhotic_associatedcells_sig)
  save(Cirrhotic_associatedcells_sig,file = "Livercirrhosis_Cirrhotic_associatedcells_sig.rdata")
  
    library(Seurat)
    library(ggplot2)
    Livercirrhosis_sc_data_downsampled$logistic_pred_labels <- factor(
      Livercirrhosis_sc_data_downsampled$logistic_pred_labels,
      levels = c(1, 0),
      labels = c("Cirrhosis-associated cells", "Other cells")
    )
    Idents(Livercirrhosis_sc_data_downsampled) <- Livercirrhosis_sc_data_downsampled@meta.data$logistic_pred_labels
    signature_genes <- rownames(Cirrhotic_associatedcells_sig)[1:10] 
    cairo_pdf(paste0("E:/project/Supervised_Phen_subpopulations_project/figures of graduation/Dot_plot_unscaled_Livercirrhosis_markers_logistic_subpopulations.pdf"), width = 10, height = 5)
    DotPlot(Livercirrhosis_sc_data_downsampled, features = signature_genes,group.by = "logistic_pred_labels",scale = FALSE) + 
      scale_color_gradient(low = "lightblue", high = "darkred") +
      theme_minimal() +
      labs(title = "Expression Levels of Signature Genes in Predicted Phenotypes",
           x = "Signature Genes",
           y = "Predicted Cell Groups") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1),
            axis.text.y = element_text(size = 10),
            axis.title = element_text(size = 14, face = "bold"),          # 坐标轴标题
            axis.line = element_line(color = "black", size = 0.8),
            panel.grid.major = element_line(color = "gray80", size = 0.5), # 主网格线
            panel.grid.minor =  element_blank()# 坐标轴线条
      )
    dev.off()
  
``````
