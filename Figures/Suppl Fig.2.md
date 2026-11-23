### Preprocessing of GSE131907 single-cell data
```r
      # ============================================================
      # GSE131907 
      # ============================================================
      
      library(Seurat)
      library(Matrix)
      library(dplyr)
      library(tibble)
      
      # ---------- 1.  ----------
      data_dir <- "/data1/liuchonghui/CTDP/data/GSE131907"
      out_dir  <- "/data1/liuchonghui/CTDP/test"

      
      anno_file <- file.path(data_dir, "GSE131907_Lung_Cancer_cell_annotation.txt")
      umi_file  <- file.path(data_dir, "GSE131907_Lung_Cancer_raw_UMI_matrix.rds")
      
      # ---------- 2.  ----------
      annotation <- read.table(
        anno_file,
        header    = TRUE,
        sep       = "\t",
        stringsAsFactors = FALSE,
        quote     = "",
        comment.char = ""
      )
      
      
      cat("注释文件维度:", dim(annotation), "\n")
      print(table(annotation$Sample_Origin))
      print(table(annotation$Cell_type))
      
      # ---------- 3.  ----------
      umi <- readRDS(umi_file)
      cat("读入对象类型:", paste(class(umi), collapse = "/"), " 维度:", dim(umi), "\n")
      
      # ---------- 4.  ----------
      mat_cells <- colnames(umi)
      head(mat_cells)
      
      anno_idx <- annotation$Index
      head(anno_idx)
      
      if (any(duplicated(anno_idx))) {
        warning("注释文件中 Index 存在重复，保留首次出现。")
        annotation <- annotation[!duplicated(anno_idx), ]
      }
      
      common_cells <- intersect(mat_cells, annotation$Index)
      cat("矩阵与注释共有细胞数:", length(common_cells), "\n")
      cat("仅在矩阵中:", length(setdiff(mat_cells, annotation$Index)), "\n")
      cat("仅在注释中:", length(setdiff(annotation$Index, mat_cells)), "\n")
      
      annotation_matched <- annotation[match(common_cells, annotation$Index), ]
      stopifnot(identical(annotation_matched$Index, colnames(umi)))
      
      metadata <- annotation_matched
      rownames(metadata) <- metadata$Index
      metadata <- metadata[, c("Barcode", "Sample", "Sample_Origin",
                               "Cell_type", "Cell_type.refined", "Cell_subtype"),
                           drop = FALSE]
      
      # ---------- 5.  ----------
      GSE131907_sc_data <- CreateSeuratObject(
        counts    = umi,
        project   = "GSE131907_Lung",
        meta.data = metadata,
        min.cells = 3,
        min.features = 200
      )
      
      cat("初始 Seurat 对象:\n")
      print(GSE131907_sc_data)
      print(table(GSE131907_sc_data$Cell_subtype, useNA = "ifany"))
      
      dim(GSE131907_sc_data)
      GSE131907_sc_data <- subset(
        GSE131907_sc_data,
        subset = grepl("^(BRONCHO|EBUS|NS|LUNG_N|LUNG_T)", Sample)
      )
      dim(GSE131907_sc_data)
      GSE131907_sc_data <- subset(
        GSE131907_sc_data,
        subset = !is.na(Cell_subtype)
      )
      dim(GSE131907_sc_data)
      GSE131907_sc_data[["percent.mt"]] <- PercentageFeatureSet(GSE131907_sc_data, pattern = "^MT-")
      GSE131907_sc_data <- subset(GSE131907_sc_data, subset = percent.mt < 20)
      dim(GSE131907_sc_data)
      GSE131907_sc_data <- NormalizeData(GSE131907_sc_data, normalization.method = "LogNormalize", scale.factor = 10000)
      GSE131907_sc_data <- ScaleData(GSE131907_sc_data)
      GSE131907_sc_data <- FindVariableFeatures(GSE131907_sc_data, selection.method = "vst", nfeatures = 2000)
      GSE131907_sc_data <- RunPCA(GSE131907_sc_data, features = VariableFeatures(object = GSE131907_sc_data))
      GSE131907_sc_data <- FindNeighbors(GSE131907_sc_data, dims = 1:10)
      GSE131907_sc_data <- FindClusters(GSE131907_sc_data)
      GSE131907_sc_data <- RunUMAP(GSE131907_sc_data, dims = 1:10)
      save(GSE131907_sc_data, file = "/data1/liuchonghui/CTDP/test/GSE131907_sc_data.rdata")
      
``````
### CTDP with cell-level permutation
```r
      load("/data1/liuchonghui/CTDP/test/GSE131907_sc_data.rdata")
      
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
        true_coefficients <- coef(lasso_model)[-1]  # 去掉截距项，保留每个 cluster 的系数
    
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
      
      meta_data <- GSE131907_sc_data@meta.data
 
      print(colnames(meta_data))
  
      df <- data.frame(
        cluster = meta_data$seurat_clusters,
        phenotype = meta_data$Sample
      )
      rownames(df) <- rownames(meta_data)

      table(df$phenotype)
      
      df$phenotype <- ifelse(
        grepl("^LUNG_N", df$phenotype), 0, 1
      )
 
      df$phenotype <- factor(df$phenotype, levels = c(0, 1))
      
      table(df$phenotype, useNA = "ifany")
  
      GSE131907_logistic_result <- lasso_logistic_regression(
        df,
        n_permutations = 1000,
        seed = 123
      )
   
      save(
        GSE131907_logistic_result,
        file = "/data1/liuchonghui/CTDP/test/GSE131907_logistic_result.rdata"
      )
``````
### CTDP with patient-level permutation
```r
      library(dplyr)
      lasso_logistic_regression_samplepermutation <- function(df, n_permutations = 1000, seed = 123) {
    
        if (!all(c("cluster", "phenotype", "sample") %in% colnames(df))) {
          stop("数据框必须包含 'cluster'、'phenotype' 和 'sample' 列。")
        }
 
        if (any(is.na(df$sample)) || any(is.na(df$phenotype))) {
          stop("sample 和 phenotype 列不能包含缺失值。")
        }
      
        set.seed(seed)
      
        x_cluster <- model.matrix(~ cluster - 1, data = df)
        x <- as.matrix(x_cluster)

        y <- df$phenotype

        cv_fit <- cv.glmnet(x, y, family = "binomial", alpha = 0)
        best_lambda <- cv_fit$lambda.min

        lasso_model <- glmnet(x, y, family = "binomial", alpha = 0, lambda = best_lambda)
        true_coefficients <- coef(lasso_model)[-1]
 
        sample_ids <- unique(df$sample)
        
        sample_phenotype <- df$phenotype[
          match(sample_ids, df$sample)
        ]
 
        sample_check <- data.frame(
          sample = df$sample,
          phenotype = df$phenotype
        )
        
        check <- sample_check %>%
          group_by(sample) %>%
          summarise(
            n_rows = n(),
            n_origin_unique = n_distinct(phenotype),
            origin_values = paste(unique(phenotype), collapse = "; "),
            .groups = "drop"
          )
        print(check, n=length(sample_ids))
        
        if (any(check$n_origin_unique != 1)) {
          problematic_samples <- check %>%
            filter(n_origin_unique != 1)
          
          print(problematic_samples)
          
          stop("检测到部分样本对应多个 phenotype，请检查上述样本的表型标签。")
        }
        
        permuted_coefficients <- matrix(
          0,
          nrow = n_permutations,
          ncol = length(true_coefficients)
        )
        
        for (i in 1:n_permutations) {

          permuted_sample_phenotype <- sample_phenotype[
            sample(seq_along(sample_phenotype))
          ]

          permuted_y <- permuted_sample_phenotype[
            match(df$sample, sample_ids)
          ]

          permuted_model <- glmnet(
            x,
            permuted_y,
            family = "binomial",
            alpha = 0,
            lambda = best_lambda
          )

          permuted_coefficients[i, ] <- coef(permuted_model)[-1]
        }

        p_values <- sapply(1:length(true_coefficients), function(j) {
          length(which(
            permuted_coefficients[, j] >= true_coefficients[j]
          )) / n_permutations
        })
        
        FDR <- p.adjust(p_values, method = "BH")

        result <- data.frame(
          Coefficient = as.vector(true_coefficients),
          P_Value = p_values,
          FDR = FDR,
          row.names = rownames(coef(lasso_model))[-1]
        )
        
        return(result)
      }
      
      load("/data1/liuchonghui/CTDP/test/GSE131907_sc_data.rdata")
      
      library(dplyr)
      library(Seurat)
      library(glmnet)
      library(SingleCellExperiment)

      meta_data <- GSE131907_sc_data@meta.data

      print(colnames(meta_data))

      df <- data.frame(
        cluster = meta_data$seurat_clusters,
        phenotype = meta_data$Sample,
        sample= meta_data$Sample
      )
      rownames(df) <- rownames(meta_data)

      table(df$phenotype)
      
      df$phenotype <- ifelse(
        grepl("^LUNG_N", df$phenotype), 0, 1
      )

      df$phenotype <- factor(df$phenotype, levels = c(0, 1))
      
      table(df$phenotype)
      
      GSE131907_logistic_result_samplepermutation <- lasso_logistic_regression_samplepermutation(
        df,
        n_permutations = 1000,
        seed = 123
      )
      
      save(
        GSE131907_logistic_result_samplepermutation,
        file = "/data1/liuchonghui/CTDP/test/GSE131907_logistic_result_samplepermutation.rdata"
      )
``````
### random scoring
```r
      # Random prediction
      set.seed(123)
      
      random_predicted_labels <- sample(
        c(0, 1),
        size = length(groundtruth_labels),
        replace = TRUE
      )
      
      table(random_predicted_labels)
      
      # Confusion matrix
      random_confusion <- confusionMatrix(
        factor(random_predicted_labels, levels = c(0, 1)),
        factor(groundtruth_labels, levels = c(0, 1)),
        positive = "1"
      )
      
      # Calculate F1-score
      random_f1_score <- random_confusion$byClass["F1"]
      
      print(random_f1_score)
``````
### plot
```r
      # ==========================================
      # Compare F1 scores of three methods
      # ==========================================
      
      library(ggplot2)
      
      # 1. Prepare data
      f1_data <- data.frame(
        Method = c(
          "CTDP",
          "CTDP with patient-level permutation",
          "Random scoring"
        ),
        F1 = c(
          0.6154477,
          0.3231908,
          0.3063283
        )
      )
      
      # Set the order of methods
      f1_data$Method <- factor(
        f1_data$Method,
        levels = c(
          "CTDP",
          "CTDP with patient-level permutation",
          "Random scoring"
        )
      )
      
      # 2. Draw bar plot
      p <- ggplot(
        f1_data,
        aes(x = Method, y = F1, fill = Method)
      ) +
        geom_bar(
          stat = "identity",
          width = 0.6,
          color = "black",
          linewidth = 0.4
        ) +
        geom_text(
          aes(label = sprintf("%.3f", F1)),
          vjust = -0.5,
          size = 5
        ) +
        scale_fill_manual(
          values = c(
            "CTDP" = "#4472C4",
            "CTDP with patient-level permutation" = "#ED7D31",
            "Random scoring" = "#A5A5A5"
          )
        ) +
        scale_y_continuous(
          limits = c(0, 0.70),
          breaks = seq(0, 0.7, by = 0.1),
          expand = expansion(mult = c(0, 0.05))
        ) +
        labs(
          x = NULL,
          y = "F1 score"
        ) +
        theme_classic(base_size = 14) +
        theme(
          legend.position = "none",
          axis.text.x = element_text(
            size = 12,
            color = "black",
            angle = 15,
            hjust = 1
          ),
          axis.text.y = element_text(
            size = 12,
            color = "black"
          ),
          axis.title.y = element_text(
            size = 14,
            color = "black"
          ),
          axis.line = element_line(
            linewidth = 0.6,
            color = "black"
          ),
          axis.ticks = element_line(
            linewidth = 0.5,
            color = "black"
          )
        )

    
      
      # 4. Save as vector PDF
      output_dir <- "/data1/liuchonghui/CTDP/figures of graduation/"
      
      if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
      }
      
      ggsave(
        filename = paste0(output_dir, "CTDP_vs_patient-level-permutation.pdf"),
        plot = p,
        device = cairo_pdf,
        width = 7,
        height = 5,
        units = "in"
      )
      
``````
