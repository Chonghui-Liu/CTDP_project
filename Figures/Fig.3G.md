## Fig3.G
### calculate F1 and generate box plot
```r
  load("/data/chliu/Supervised_Phen_subpopulations_project/test/PBMC_Travaglini_SC_sample3clusters_list.rdata")
  
  
  load("/data/chliu/Supervised_Phen_subpopulations_project/test/logistic_result_list_sample3clusters.rdata")
  library(dplyr)
  library(caret)
  logistic_f1_vector <- numeric(length(logistic_result_list_sample3clusters))
  for (i in seq_along(logistic_result_list_sample3clusters)) {
    logistic_result <- logistic_result_list_sample3clusters[[i]]
    seurat_obj <- PBMC_Travaglini_SC_sample3clusters_list[[i]]
    significant_clusters <- rownames(logistic_result[logistic_result$FDR < 0.05, ])
    pred_labels <- as.data.frame(
      ifelse(seurat_obj@meta.data$seurat_clusters %in% as.numeric(gsub("cluster", "", significant_clusters)), 1, 0)
    )
    colnames(pred_labels) <- "logistic_pred_labels"
    rownames(pred_labels) <- rownames(seurat_obj@meta.data)
    seurat_obj <- AddMetaData(seurat_obj, metadata = pred_labels)
    groundtruth_labels <- seurat_obj@meta.data$groundtruth
    predicted_labels <- seurat_obj@meta.data$logistic_pred_labels
    confusion <- confusionMatrix(factor(predicted_labels), factor(groundtruth_labels))
    f1_score <- confusion$byClass["F1"]
    logistic_f1_vector[i] <- f1_score
  }
  print(logistic_f1_vector)
  
  
  load("/data/chliu/Supervised_Phen_subpopulations_project/test/scDist_result_list_sample3clusters.rdata")
  scDist_f1_vector <- numeric(length(scDist_result_list_sample3clusters))
  for (i in seq_along(scDist_result_list_sample3clusters)) {
    scDist_result <- scDist_result_list_sample3clusters[[i]]
    seurat_obj <- PBMC_Travaglini_SC_sample3clusters_list[[i]]
    significant_clusters <- rownames(scDist_result$results[scDist_result$results$p.val < 0.05, ])
    print(significant_clusters)
    pred_labels <- as.data.frame(
      ifelse(seurat_obj@meta.data$seurat_clusters %in% as.numeric(significant_clusters), 1, 0)
    )
    colnames(pred_labels) <- "scDist_pred_labels"
    rownames(pred_labels) <- rownames(seurat_obj@meta.data)
    print(table(pred_labels))
    seurat_obj <- AddMetaData(seurat_obj, metadata = pred_labels)
    groundtruth_labels <- seurat_obj@meta.data$groundtruth
    predicted_labels <- seurat_obj@meta.data$scDist_pred_labels
    confusion <- confusionMatrix(factor(predicted_labels), factor(groundtruth_labels))
    f1_score <- confusion$byClass["F1"]
    scDist_f1_vector[i] <- f1_score
  }
  print(scDist_f1_vector)
  
  
  load("/data/chliu/Supervised_Phen_subpopulations_project/test/da_cells_list_sample3clusters.rdata")
  DAseq_f1_vector <- lapply(seq_along(PBMC_Travaglini_SC_sample3clusters_list), function(i) {
    # 获取当前的 Seurat 对象和 DA 细胞检测结果
    sc_data <- PBMC_Travaglini_SC_sample3clusters_list[[i]]
    da_cells <- da_cells_list_sample3clusters[[i]]  # 对应的 DA 细胞结果
    
    # 初始化预测标签
    pred_labels <- rep(0, nrow(sc_data@meta.data))
    pred_labels[da_cells$da.down] <- 1
    pred_labels <- as.data.frame(pred_labels)
    colnames(pred_labels) <- "DAseq_pred_labels"
    rownames(pred_labels) <- rownames(sc_data@meta.data)
    print(table(pred_labels))
    
    # 添加预测标签到元数据
    sc_data <- AddMetaData(sc_data, metadata = pred_labels)
    
    # 提取真实标签和预测标签
    groundtruth_labels <- sc_data@meta.data$groundtruth
    predicted_labels <- sc_data@meta.data$DAseq_pred_labels
    
    # 计算混淆矩阵和 F1 得分
    confusion <- confusionMatrix(factor(predicted_labels), factor(groundtruth_labels))
    f1_score <- confusion$byClass["F1"]
    
    return(f1_score)
  })
  DAseq_f1_vector <- unlist(DAseq_f1_vector)
  print(DAseq_f1_vector)
  
  
  PENCIL_f1_vector <- numeric(length(PBMC_Travaglini_SC_sample3clusters_list))
  library(Seurat)
  PENCIL_f1_vector <- sapply(seq_along(PBMC_Travaglini_SC_sample3clusters_list), function(i) {
    sc_data <- PBMC_Travaglini_SC_sample3clusters_list[[i]]
    
    # 导入Python生成的预测标签
    pred_labels <- read.csv(paste0("/data/chliu/Supervised_Phen_subpopulations_project/test/PBMC_Travaglini_SC_3clusterslabeled_", i-1, "_PENCILpred_labels.csv"), header = FALSE)
    colnames(pred_labels) <- "PENCIL_pred_labels"
    rownames(pred_labels) <- rownames(sc_data@meta.data)
    
    # 将预测标签添加到Seurat对象
    sc_data <- AddMetaData(sc_data, metadata = pred_labels)
    
    # 计算F1分数
    groundtruth_labels <- sc_data@meta.data$groundtruth
    predicted_labels <- sc_data@meta.data$PENCIL_pred_labels
    predicted_numeric <- ifelse(predicted_labels == "1.0", 1, 0)
    confusion <- confusionMatrix(factor(predicted_numeric), factor(groundtruth_labels))
    f1_score <- confusion$byClass["F1"]
    
    return(f1_score)
  })
  print(PENCIL_f1_vector)
  
  
  library(ggplot2)
  f1_data <- data.frame(
    F1_Score = c(logistic_f1_vector, scDist_f1_vector, DAseq_f1_vector, PENCIL_f1_vector),
    Method = factor(rep(c("Logistic", "scDist", "DAseq", "PENCIL"), each = length(logistic_f1_vector)))
  )
  cairo_pdf("/data/chliu/Supervised_Phen_subpopulations_project/figures of graduation/methods_comparision_sample3clusters.pdf")
  ggplot(f1_data, aes(x = Method, y = F1_Score, fill = Method)) +
    geom_boxplot() +
    labs(title = "Comparison of F1 Scores across Different Methods",
         x = "Method",
         y = "F1 Score") +
    theme_minimal() +
    theme(
      legend.position = "none",
      panel.border = element_rect(color = "black", fill = NA, size = 1),# 添加四周的黑色边框
      panel.grid = element_blank()
    ) +
    scale_fill_brewer(palette = "Set3")
  dev.off()


``````
