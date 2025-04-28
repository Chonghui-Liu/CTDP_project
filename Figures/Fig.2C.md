## Fig2.C
### Cell clusters 0-9 were labeled with mock labels, respectively
```r
  assign_labels <- function(seurat_obj, clusters, ratio) {
    meta_data <- seurat_obj@meta.data
    target_cells <- which(meta_data$seurat_clusters %in% clusters)
    
    meta_data$source_phen <- NA
    meta_data$groundtruth <- 0
    
    meta_data$groundtruth[target_cells] <- 1
    meta_data$source_phen[target_cells] <- ifelse(runif(length(target_cells)) < (1 - ratio), 1, 0)
    
    non_target_cells <- setdiff(1:nrow(meta_data), target_cells)
    meta_data$source_phen[non_target_cells] <- sample(c(0, 1), length(non_target_cells), replace = TRUE)
    
    seurat_obj@meta.data <- meta_data
    return(seurat_obj)
  }
  
  load("/data/chliu/Supervised_Phen_subpopulations_project/test/PBMC_Travaglini_SC.rdata")
  PBMC_Travaglini_SC_cluster0to9labeled_list <- list()
  for (cluster in 0:9) {
    labeled_obj <- assign_labels(PBMC_Travaglini_SC, clusters = c(cluster), ratio = 0.1)
    name <- paste0("PBMC_Travaglini_SC_cluster", cluster, "labeled")
    PBMC_Travaglini_SC_cluster0to9labeled_list[[name]] <- labeled_obj
  }
  save(PBMC_Travaglini_SC_cluster0to9labeled_list,file = "/data/chliu/Supervised_Phen_subpopulations_project/test/PBMC_Travaglini_SC_cluster0to9labeled_list.rdata")  

``````
### The CTDP algorithm is executed
```r
  load("/data/chliu/Supervised_Phen_subpopulations_project/test/PBMC_Travaglini_SC_cluster0to9labeled_list.rdata")
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
  logistic_result_list <- lapply(PBMC_Travaglini_SC_cluster0to9labeled_list, function(seurat_obj) {
    df <- data.frame(
      cluster = seurat_obj@meta.data$seurat_clusters,
      phenotype = seurat_obj@meta.data$source_phen
    )
    lasso_logistic_regression(df, n_permutations = 1000, seed = 123)
  })
  save(logistic_result_list,file = "/data/chliu/Supervised_Phen_subpopulations_project/test/logistic_result_list.rdata")
      
  
``````
### Generate an UMAP plot
```r
  load("/data/chliu/Supervised_Phen_subpopulations_project/test/PBMC_Travaglini_SC_cluster0to9labeled_list.rdata")
  i = 2
  print(names(PBMC_Travaglini_SC_cluster0to9labeled_list)[i])
  seurat_obj <- PBMC_Travaglini_SC_cluster0to9labeled_list[[i]]
  
  
  load("/data/chliu/Supervised_Phen_subpopulations_project/test/logistic_result_list.rdata")
  logistic_result <- logistic_result_list[[i]]
  significant_clusters <- rownames(logistic_result[logistic_result$FDR < 0.05, ])
  print(significant_clusters)
  pred_labels <- as.data.frame(
    ifelse(seurat_obj@meta.data$seurat_clusters %in% as.numeric(gsub("cluster", "", significant_clusters)), 1, 0)
  )
  table(pred_labels)
  colnames(pred_labels) <- "logistic_pred_labels"
  rownames(pred_labels) <- rownames(seurat_obj@meta.data)
  seurat_obj <- AddMetaData(seurat_obj, metadata = pred_labels)
  seurat_obj@meta.data[1:5,16:ncol(seurat_obj@meta.data)]

  cairo_pdf(paste0("/data/chliu/Supervised_Phen_subpopulations_project/figures of graduation/",names(PBMC_Travaglini_SC_cluster0to9labeled_list)[i],"_logistic_pred_UMAP.pdf"))
  DimPlot(seurat_obj, group.by = "logistic_pred_labels", reduction = "umap", label = F, repel = TRUE,pt.size = 0.5) + NoLegend()+theme(aspect.ratio=1, legend.position = "right")
  dev.off()
  
``````