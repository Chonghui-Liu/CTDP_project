## Fig4.D
```r
  load("/data/chliu/Supervised_Phen_subpopulations_project/test/GSE120575_logistic_result.rdata")
  load("/data/chliu/Supervised_Phen_subpopulations_project/test/GSE120575_sc_data_downsampled.rdata")
  logistic_result <- GSE120575_logistic_result
  significant_clusters <- rownames(logistic_result[logistic_result$FDR < 0.05, ])
  print(significant_clusters)
  pred_labels <- as.data.frame(
    ifelse(GSE120575_sc_data_downsampled@meta.data$seurat_clusters %in% as.numeric(gsub("cluster", "", significant_clusters)), 1, 0)
  )
  table(pred_labels)
  colnames(pred_labels) <- "logistic_pred_labels"
  rownames(pred_labels) <- rownames(GSE120575_sc_data_downsampled@meta.data)
  GSE120575_sc_data_downsampled <- AddMetaData(GSE120575_sc_data_downsampled, metadata = pred_labels)
  
  Responder_associatedcells_sig <- Seurat::FindMarkers(GSE120575_sc_data_downsampled, ident.1 = 1, group.by = "logistic_pred_labels", logfc.threshold=0, min.pct=0,random.seed=123)
  Responder_associatedcells_sig$names <- rownames(Responder_associatedcells_sig)
  save(Responder_associatedcells_sig,file = "/data/chliu/Supervised_Phen_subpopulations_project/test/GSE120575_Responder_associatedcells_sig.rdata")
  
  
  library(Seurat)
  library(ggplot2)
  Idents(GSE120575_sc_data_downsampled) <- GSE120575_sc_data_downsampled@meta.data$response
  signature_genes <- rownames(Responder_associatedcells_sig)[1:10] 
  cairo_pdf(paste0("/data/chliu/Supervised_Phen_subpopulations_project/figures of graduation/Dot_plot_GSE120575_markers_logistic_subpopulations.pdf"), width = 10, height = 5)
  DotPlot(GSE120575_sc_data_downsampled, features = signature_genes,group.by = "response") + 
    scale_color_gradient2(low = "lightblue", high = "darkred", mid = "white", midpoint = 0) +
    theme_minimal() +
    labs(title = "Expression Levels of Signature Genes in Predicted Phenotypes",
         x = "Signature Genes",
         y = "Response Group") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          axis.text.y = element_text(size = 10),
          axis.title = element_text(size = 14, face = "bold"),         
          axis.line = element_line(color = "black", size = 0.8),
          panel.grid.major = element_line(color = "gray80", size = 0.5), 
          panel.grid.minor =  element_blank()
          )
  dev.off()  
  
``````
