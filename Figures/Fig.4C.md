## Fig4.C
### Perform the CTDP method on the single-cell data of melanoma cases
```r
  load("/data/chliu/Supervised_Phen_subpopulations_project/test/GSE120575_sc_data_downsampled.rdata")
  library(dplyr)
  library(glmnet)
  lasso_logistic_regression <- function(df, n_permutations = 1000, seed = 123) {
    if (!all(c("cluster", "phenotype") %in% colnames(df))) {
      stop("data_frame must contain 'cluster' and 'phenotype' column。")
    }
    
    set.seed(seed)
    
    x_cluster <- model.matrix(~ cluster - 1, data = df)
    x <- as.matrix(x_cluster)
    
    y <- df$phenotype
    
    cv_fit <- cv.glmnet(x, y, family = "binomial", alpha = 0)
    best_lambda <- cv_fit$lambda.min
    
    lasso_model <- glmnet(x, y, family = "binomial", alpha = 0, lambda = best_lambda)
    true_coefficients <- coef(lasso_model)[-1]  
    
    permuted_coefficients <- matrix(0, nrow = n_permutations, ncol = length(true_coefficients))
    for (i in 1:n_permutations) {
      permuted_y <- sample(y)  
      permuted_model <- glmnet(x, permuted_y, family = "binomial", alpha = 0, lambda = best_lambda)
      permuted_coefficients[i, ] <- coef(permuted_model)[-1]  
    }
    
    p_values <- sapply(1:length(true_coefficients), function(j) {
      length(which(permuted_coefficients[, j] >= true_coefficients[j]))/n_permutations
    })
    
    FDR <- p.adjust(p_values,method = "BH")
    result <- data.frame(
      Coefficient = as.vector(true_coefficients),
      P_Value = p_values,
      FDR=FDR,
      row.names = rownames(coef(lasso_model))[-1]  
    )
    
    return(result)
  }
  df <- data.frame(
    cluster = GSE120575_sc_data_downsampled@meta.data$seurat_clusters,
    phenotype = GSE120575_sc_data_downsampled@meta.data$response
  )
  GSE120575_logistic_result <- lasso_logistic_regression(df, n_permutations = 1000, seed = 123)
  save(GSE120575_logistic_result,file = "/data/chliu/Supervised_Phen_subpopulations_project/test/GSE120575_logistic_result.rdata")
  
``````
### Plot Fig.4C
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
  GSE120575_sc_data_downsampled@meta.data[1:5,]
  
  cairo_pdf(paste0("/data/chliu/Supervised_Phen_subpopulations_project/figures of graduation/GSE120575_sc_data_downsampled_logistic_pred_UMAP.pdf"))
  DimPlot(GSE120575_sc_data_downsampled, group.by = "logistic_pred_labels", reduction = "umap", label = F, repel = TRUE,pt.size = 0.5) + NoLegend()+theme(aspect.ratio=1, legend.position = "right")
  dev.off()  
  
``````