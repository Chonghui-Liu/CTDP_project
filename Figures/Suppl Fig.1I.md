### Draw a box plot comparing several methods based on GSE182227 data
```r
# ============================================================
      # 1. load data
      # ============================================================
      
      library(dplyr)
      library(caret)
      library(Seurat)
      library(ggplot2)
      
      load("/data1/liuchonghui/CTDP/test/GSE182227_sc_data.rdata")
      
      meta_data <- GSE182227_sc_data@meta.data
      
      groundtruth_labels <- meta_data$Cancer_assignment
      
      table(groundtruth_labels)
      
      groundtruth_labels <- factor(
        ifelse(groundtruth_labels %in% c("Cancer"), 1, 0),
        levels = c(1, 0)
      )
      
      table(groundtruth_labels)
      
      
      # ============================================================
      # 2. Logistic 
      # ============================================================
      
      load("/data1/liuchonghui/CTDP/test/GSE182227_logistic_result.rdata")
      
      logistic_result <- GSE182227_logistic_result
      seurat_obj <- GSE182227_sc_data
      
      significant_clusters <- rownames(
        logistic_result[logistic_result$FDR < 0.05, ]
      )
      significant_clusters
      
      pred_labels <- as.data.frame(
        ifelse(
          meta_data$seurat_clusters %in% as.numeric(gsub("cluster", "", significant_clusters)),
          1,
          0
        )
      )
      
      colnames(pred_labels) <- "logistic_pred_labels"
      rownames(pred_labels) <- rownames(meta_data)
      
      table(pred_labels)
      head(pred_labels)
      
      predicted_labels <- pred_labels$logistic_pred_labels
      
      confusion <- confusionMatrix(
        factor(predicted_labels, levels = c(0, 1)),
        factor(groundtruth_labels, levels = c(0, 1)),
        positive = "1"
      )
      
      logistic_f1 <- confusion$byClass["F1"]
      
      print(logistic_f1)
      
      
      # ============================================================
      # 3. DAseq 
      # ============================================================
      
      load("/data1/liuchonghui/CTDP/test/GSE182227_DAseq_result.rdata")
      
      da_cells <- GSE182227_DAseq_result
     
      pred_labels <- rep(0, nrow(meta_data))
      pred_labels[da_cells$da.down] <- 1
      
      pred_labels <- as.data.frame(pred_labels)
      colnames(pred_labels) <- "DAseq_pred_labels"
      rownames(pred_labels) <- rownames(meta_data)
      
      print(table(pred_labels))
      
      predicted_labels <- pred_labels$DAseq_pred_labels
  
      valid_cells <- !is.na(groundtruth_labels)
      
      DAseq_groundtruth_labels <- groundtruth_labels[valid_cells]
      DAseq_predicted_labels <- predicted_labels[valid_cells]

      confusion <- confusionMatrix(
        factor(DAseq_predicted_labels, levels = c(0, 1)),
        factor(DAseq_groundtruth_labels, levels = c(0, 1)),
        positive = "1"
      )
      
      DAseq_f1 <- confusion$byClass["F1"]
      
      print(DAseq_f1)
      
      
      # ============================================================
      # 4. scDist 
      # ============================================================
      
      load("/data1/liuchonghui/CTDP/test/GSE182227_scDist_result.rdata")
      
      significant_clusters <- rownames(
        GSE182227_scDist_result$results[
          GSE182227_scDist_result$results$p.val < 0.05,
        ]
      )
      
      print(significant_clusters)
      
      pred_labels <- as.data.frame(
        ifelse(
          meta_data$seurat_clusters %in%
            as.numeric(significant_clusters),
          1,
          0
        )
      )
      
      colnames(pred_labels) <- "scDist_pred_labels"
      rownames(pred_labels) <- rownames(meta_data)
      
      print(table(pred_labels))
      
      GSE182227_sc_data <- AddMetaData(
        GSE182227_sc_data,
        metadata = pred_labels
      )
      
      predicted_labels <- GSE182227_sc_data@meta.data$scDist_pred_labels
      
      confusion <- confusionMatrix(
        factor(predicted_labels, levels = c(0, 1)),
        factor(groundtruth_labels, levels = c(0, 1)),
        positive = "1"
      )
      
      print(confusion)
      
      scDist_f1 <- confusion$byClass["F1"]
      
      print(scDist_f1)
      
      
      # ============================================================
      # 5. PENCIL 
      # ===========================================================
      
      pred_labels <- read.csv(
        "/data1/liuchonghui/CTDP/test/GSE182227_sc_data_scaled_PENCILpred_labels.csv",
        header = FALSE
      )
      
      colnames(pred_labels) <- "PENCIL_pred_labels"
      
      head(pred_labels)
      table(pred_labels)
      
      rownames(pred_labels) <- rownames(meta_data)
      
      predicted_labels <- pred_labels$PENCIL_pred_labels
      
      stopifnot(length(predicted_labels) == length(groundtruth_labels))
      
      table(groundtruth_labels)
      table(predicted_labels)
  
      predicted_numeric <- ifelse(predicted_labels == "1", 1, 0)
      
      table(predicted_numeric)

      confusion <- confusionMatrix(
        factor(predicted_numeric, levels = c(0, 1)),
        factor(groundtruth_labels, levels = c(0, 1)),
        positive = "1"
      )

      PENCIL_f1 <- confusion$byClass["F1"]
      
      print(PENCIL_f1)
      

      # ============================================================
      # 6. Chi-square 
      # ============================================================
      load("/data1/liuchonghui/CTDP/test/GSE182227_chi_square_result.rdata")
      
      chi_square_result <- GSE182227_chi_square_result
      
      significant_clusters <- rownames(
        chi_square_result[chi_square_result$FDR < 0.05, ]
      )
      
      print(significant_clusters)
      
      pred_labels <- as.data.frame(
        ifelse(
          meta_data$seurat_clusters %in%
            as.numeric(gsub("cluster", "", significant_clusters)),
          1,
          0
        )
      )
      
      colnames(pred_labels) <- "chi_square_pred_labels"
      rownames(pred_labels) <- rownames(meta_data)
      
      print(table(pred_labels))
      
      predicted_labels <- pred_labels$chi_square_pred_labels
      
      confusion <- confusionMatrix(
        factor(predicted_labels, levels = c(0, 1)),
        factor(groundtruth_labels, levels = c(0, 1)),
        positive = "1"
      )
      
      chi_square_f1 <- confusion$byClass["F1"]
      
      print(chi_square_f1)
      
      
      # ============================================================
      # 7. Pearson 
      # ============================================================
      load("/data1/liuchonghui/CTDP/test/GSE182227_pearson_result.rdata")
      
      pearson_result <- GSE182227_pearson_result
      
      significant_clusters <- rownames(
        pearson_result[pearson_result$FDR < 0.05, ]
      )
      
      print(significant_clusters)
      
      pred_labels <- as.data.frame(
        ifelse(
          meta_data$seurat_clusters %in%
            as.numeric(gsub("cluster", "", significant_clusters)),
          1,
          0
        )
      )
      
      colnames(pred_labels) <- "pearson_pred_labels"
      rownames(pred_labels) <- rownames(meta_data)
      
      print(table(pred_labels))
      
      predicted_labels <- pred_labels$pearson_pred_labels
      
      confusion <- confusionMatrix(
        factor(predicted_labels, levels = c(0, 1)),
        factor(groundtruth_labels, levels = c(0, 1)),
        positive = "1"
      )
      
      pearson_f1 <- confusion$byClass["F1"]
      
      print(pearson_f1)
      
      
      # ============================================================
      # 8. summary 
      # ============================================================
      
      f1_data <- data.frame(
        Method = factor(
          c("Logistic", "Chi-square", "Pearson",
            "scDist", "DAseq", "PENCIL"),
          levels = c("Logistic", "Chi-square", "Pearson",
                     "scDist", "DAseq", "PENCIL")
        ),
        F1_Score = c(
          as.numeric(logistic_f1),
          as.numeric(chi_square_f1),
          as.numeric(pearson_f1),
          as.numeric(scDist_f1),
          as.numeric(DAseq_f1),
          as.numeric(PENCIL_f1)
        )
      )
      
      print(f1_data)
      
      
      # ============================================================
      # ============================================================
      # 9. plot
      # ============================================================
      
      p <- ggplot(
        f1_data,
        aes(x = F1_Score, y = Method, color = Method)
      ) +
        geom_segment(
          aes(
            x = 0,
            xend = F1_Score,
            y = Method,
            yend = Method
          ),
          color = "grey75",
          linewidth = 0.8
        ) +
        geom_point(size = 4) +
        geom_text(
          aes(label = sprintf("%.3f", F1_Score)),
          hjust = -0.5,
          size = 4,
          show.legend = FALSE
        ) +
        scale_x_continuous(
          limits = c(0, 1.08),
          breaks = seq(0, 1, 0.2),
          expand = c(0, 0)
        ) +
        scale_color_brewer(palette = "Dark2") +
        labs(
          x = "F1-score",
          y = NULL
        ) +
        theme_classic(base_size = 14) +
        theme(
          legend.position = "none",
          axis.text = element_text(
            color = "black",
            size = 12
          ),
          axis.title.x = element_text(
            color = "black",
            size = 14
          ),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.border = element_rect(
            color = "black",
            fill = NA,
            linewidth = 0.8
          )
        )
      

      ggsave(
        filename = "/data1/liuchonghui/CTDP/figures of graduation/GSE182227_six_methods_comparison_F1.pdf",
        plot = p,
        width = 6,
        height = 5,
        device = cairo_pdf
      )
``````
