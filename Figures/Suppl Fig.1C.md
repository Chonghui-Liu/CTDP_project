### The CTDP algorithm is executed
```r
      load("/data1/liuchonghui/CTDP/test/GSE182227_sc_data.rdata")
      
      library(dplyr)
      library(Seurat)
      library(glmnet)
      library(SingleCellExperiment)
      
      lasso_logistic_regression <- function(df, n_permutations = 1000, seed = 123) {
        if (!all(c("cluster", "phenotype") %in% colnames(df))) {
          stop("数据框必须包含 'cluster' 和 'phenotype' 列。")
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
          row.names = rownames(coef(lasso_model))[-1]  # 设置行为 cluster 名称
        )
        
        return(result)
      }

      meta_data <- GSE182227_sc_data@meta.data
      

      print(colnames(meta_data))
      

      df <- data.frame(
        cluster = meta_data$seurat_clusters,
        phenotype = meta_data$Sequenced_sample
      )
      rownames(df) <- rownames(meta_data)
      

      table(df$phenotype)
      df$phenotype <- factor(
        ifelse(df$phenotype %in% c("OP33", "OP34", "OP35BOT"), 1, 0),
        levels = c(0, 1)
      )
      table(df$phenotype)
      
      GSE182227_logistic_result <- lasso_logistic_regression(
        df,
        n_permutations = 1000,
        seed = 123
      )
      
      save(
        GSE182227_logistic_result,
        file = "/data1/liuchonghui/CTDP/test/GSE182227_logistic_result.rdata"
      )
      
``````
### Generate an UMAP plot
```r
      load("/data1/liuchonghui/CTDP/test/GSE182227_sc_data.rdata")
      seurat_obj <- GSE182227_sc_data
      

      load("/data1/liuchonghui/CTDP/test/GSE182227_logistic_result.rdata")
      
      logistic_result <- GSE182227_logistic_result
      
      significant_clusters <- rownames(
        logistic_result[logistic_result$FDR < 0.05, ]
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
      
      colnames(pred_labels) <- "logistic_pred_labels"
      rownames(pred_labels) <- rownames(seurat_obj@meta.data)
      
      seurat_obj <- AddMetaData(seurat_obj, metadata = pred_labels)
      seurat_obj@meta.data[1:5, ncol(seurat_obj@meta.data)-5:ncol(seurat_obj@meta.data)]

      library(Seurat)
      library(ggpubr)
      library(ggplot2)
      cairo_pdf(
        "/data1/liuchonghui/CTDP/figures of graduation/GSE182227_logistic_pred_UMAP.pdf"
      )
      DimPlot(
        seurat_obj,
        group.by = "logistic_pred_labels",
        reduction = "umap",
        label = F,
        repel = TRUE,
        pt.size = 0.5
      ) + NoLegend() +
        theme(aspect.ratio = 1, legend.position = "right")
      dev.off()

``````
