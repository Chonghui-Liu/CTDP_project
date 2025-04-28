## Fig5.D
```r
load("/data/chliu/Supervised_Phen_subpopulations_project/test/COVID19_logistic_result.rdata")
  load("/data/chliu/Supervised_Phen_subpopulations_project/test/COVID19_sc_data_downsampled.rdata")
  logistic_result <- COVID19_logistic_result
  significant_clusters <- rownames(logistic_result[logistic_result$FDR < 0.05, ])
  print(significant_clusters)
  pred_labels <- as.data.frame(
    ifelse(COVID19_sc_data_downsampled@meta.data$celltype %in% as.character(gsub("cluster", "", significant_clusters)), 1, 0)
  )
  table(pred_labels)
  colnames(pred_labels) <- "logistic_pred_labels"
  rownames(pred_labels) <- rownames(COVID19_sc_data_downsampled@meta.data)
  COVID19_sc_data_downsampled <- AddMetaData(COVID19_sc_data_downsampled, metadata = pred_labels)
  
  options(future.globals.maxSize = 5 * 1024^3)
  library(future)
  plan("multicore", workers = 50)
  critical_associatedcells_sig <- Seurat::FindMarkers(COVID19_sc_data_downsampled, ident.1 = 1, group.by = "logistic_pred_labels", logfc.threshold=0, min.pct=0,random.seed=123)
  critical_associatedcells_sig$names <- rownames(critical_associatedcells_sig)
  save(critical_associatedcells_sig,file = "/data/chliu/Supervised_Phen_subpopulations_project/test/COVID19_critical_associatedcells_sig.rdata")
  
  
  library(Seurat)
  library(ggplot2)
  Idents(COVID19_sc_data_downsampled) <- COVID19_sc_data_downsampled@meta.data$severity
  signature_genes <- rownames(critical_associatedcells_sig)[1:10] 
  cairo_pdf(paste0("/data/chliu/Supervised_Phen_subpopulations_project/figures of graduation/Dot_plot_COVID19_markers_logistic_subpopulations.pdf"), width = 10, height = 5)
  DotPlot(COVID19_sc_data_downsampled, features = signature_genes,group.by = "severity") + 
    scale_color_gradient2(low = "lightblue", high = "darkred", mid = "white", midpoint = 0) +
    theme_minimal() +
    labs(title = "Expression Levels of Signature Genes in Predicted Phenotypes",
         x = "Signature Genes",
         y = "severity Group") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          axis.text.y = element_text(size = 10),
          axis.title = element_text(size = 14, face = "bold"),          
          axis.line = element_line(color = "black", size = 0.8),
          panel.grid.major = element_line(color = "gray80", size = 0.5), 
          panel.grid.minor =  element_blank()
    )
  dev.off()  
  
``````
