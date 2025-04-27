## Fig6.C
### Perform the CTDP method on the single-cell data of liver cirrhosis cases
```r
  load("Livercirrhosis_sc_data_downsampled.rdata")
  library(dplyr)
  library(glmnet)
  CTDP <- function(df, n_permutations = 1000, seed = 123) {
    if (!all(c("cluster", "phenotype") %in% colnames(df))) {
      stop("data_frame must contain 'cluster' 和 'phenotype' column。")
    }
    
    set.seed(seed)
    
    x_cluster <- model.matrix(~ cluster - 1, data = df)
    x <- as.matrix(x_cluster)
    
    y<- factor(df$phenotype, levels = c("Uninjured","Cirrhotic"))
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
    cluster = Livercirrhosis_sc_data_downsampled@meta.data$annotation_lineage,
    phenotype = Livercirrhosis_sc_data_downsampled@meta.data$condition
  )
  head(df)
  levels(factor(df$phenotype))
  Livercirrhosis_logistic_result <- CTDP(df, n_permutations = 1000, seed = 123)
  Livercirrhosis_logistic_result
  save(Livercirrhosis_logistic_result,file = "Livercirrhosis_logistic_result.rdata")

``````
### Plot Fig.6C
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
  
  
  cairo_pdf(paste0("Livercirrhosis_sc_data_downsampled_logistic_pred_UMAP.pdf"))
  DimPlot(Livercirrhosis_sc_data_downsampled, group.by = "logistic_pred_labels", reduction = "umap", label = F, repel = TRUE,pt.size = 0.5) + NoLegend()+theme(aspect.ratio=1, legend.position = "right")
  dev.off()
   
``````
