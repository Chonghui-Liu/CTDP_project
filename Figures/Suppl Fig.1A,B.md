## Supplementary Fig1.A,B
### Preprocessing of GSE182227 single-cell data
```r
      library(Matrix)
      library(Seurat)
      
      data_dir <- "/data1/liuchonghui/CTDP/data/GSE182227"
      output_dir <- "/data1/liuchonghui/CTDP/test"
      
      expr_matrix <- readMM(
        file.path(data_dir, "GSE182227_OPSCC.mtx")
      )
      dim(expr_matrix)
      
      gene_names <- readLines(
        file.path(data_dir, "GSE182227_genes.txt")
      )
      
      barcodes <- readLines(
        file.path(data_dir, "GSE182227_barcodes.txt")
      )
      
      print(dim(expr_matrix))
      print(length(gene_names))
      print(length(barcodes))
      
      stopifnot(
        nrow(expr_matrix) == length(gene_names),
        ncol(expr_matrix) == length(barcodes)
      )
      
      rownames(expr_matrix) <- make.unique(gene_names)
      colnames(expr_matrix) <- make.unique(barcodes)

      GSE182227_sc_data <- CreateSeuratObject(
        counts = expr_matrix,
        project = "GSE182227_OPSCC",
        min.cells = 3,
        min.features = 200
      )
      
      cell_annotation <- read.delim(
        file.path(data_dir, "GSE182227_cell_annotation.tsv"),
        header = TRUE,
        sep = "\t",
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
      
      print(colnames(cell_annotation))
      print(head(cell_annotation))

      barcode_col <- "Cell_barcode"
      
      stopifnot(barcode_col %in% colnames(cell_annotation))
      
      stopifnot(!anyDuplicated(cell_annotation[[barcode_col]]))
      
      rownames(cell_annotation) <- cell_annotation[[barcode_col]]
      
      cell_annotation[[barcode_col]] <- NULL
      
      seurat_cells <- colnames(GSE182227_sc_data)
      annotation_cells <- rownames(cell_annotation)
      
      cat("表达矩阵细胞数：", length(seurat_cells), "\n")
      cat("注释文件细胞数：", length(annotation_cells), "\n")
      cat("成功匹配细胞数：",
          sum(seurat_cells %in% annotation_cells), "\n")
      cat("未匹配注释的细胞数：",
          sum(!seurat_cells %in% annotation_cells), "\n")
      cat("注释文件中未匹配表达矩阵的细胞数：",
          sum(!annotation_cells %in% seurat_cells), "\n")
      all(seurat_cells==annotation_cells)
      
      cell_annotation_matched <- cell_annotation[
        match(seurat_cells, annotation_cells),
        ,
        drop = FALSE
      ]
      
      rownames(cell_annotation_matched) <- seurat_cells
      
      GSE182227_sc_data <- AddMetaData(
        object = GSE182227_sc_data,
        metadata = cell_annotation_matched
      )
      
      print(dim(GSE182227_sc_data))
      print(colnames(GSE182227_sc_data@meta.data))
      print(head(GSE182227_sc_data@meta.data))
      
      
      selected_samples <- c(
        "OP33",
        "OP33norm",
        "OP34",
        "OP34norm",
        "OP35BOT",
        "OP35norm"
      )
      
      stopifnot("orig.ident" %in% colnames(GSE182227_sc_data@meta.data))
      
      print(table(GSE182227_sc_data$orig.ident))
      
      GSE182227_sc_data_subset <- subset(
        GSE182227_sc_data,
        subset = orig.ident %in% selected_samples
      )
      
      print(GSE182227_sc_data_subset)
      print(table(GSE182227_sc_data_subset$orig.ident))
      
      missing_samples <- setdiff(
        selected_samples,
        unique(as.character(GSE182227_sc_data_subset$orig.ident))
      )
      
      if (length(missing_samples) > 0) {
        warning(
          "以下样本未在提取结果中找到：",
          paste(missing_samples, collapse = ", ")
        )
      }
      
      GSE182227_sc_data = GSE182227_sc_data_subset
      dim(GSE182227_sc_data)
      GSE182227_sc_data[["percent.mt"]] <- PercentageFeatureSet(GSE182227_sc_data, pattern = "^MT-")
      GSE182227_sc_data <- subset(GSE182227_sc_data, subset = percent.mt < 20)
      dim(GSE182227_sc_data)
      GSE182227_sc_data <- NormalizeData(GSE182227_sc_data, normalization.method = "LogNormalize", scale.factor = 10000)
      GSE182227_sc_data <- ScaleData(GSE182227_sc_data)
      GSE182227_sc_data <- FindVariableFeatures(GSE182227_sc_data, selection.method = "vst", nfeatures = 2000)
      GSE182227_sc_data <- RunPCA(GSE182227_sc_data, features = VariableFeatures(object = GSE182227_sc_data))
      GSE182227_sc_data <- FindNeighbors(GSE182227_sc_data, dims = 1:10)
      GSE182227_sc_data <- FindClusters(GSE182227_sc_data)
      GSE182227_sc_data <- RunUMAP(GSE182227_sc_data, dims = 1:10)
      save(GSE182227_sc_data,file="/data1/liuchonghui/CTDP/test/GSE182227_sc_data.rdata") 
``````
### Generate an UMAP plot
```r
      load("/data1/liuchonghui/CTDP/test/GSE182227_sc_data.rdata")
      seurat_obj <- GSE182227_sc_data

      library(Seurat)
      library(ggpubr)
      library(ggplot2)
      
      cairo_pdf(
        "/data1/liuchonghui/CTDP/figures of graduation/GSE182227_UMAP.pdf"
      )
      DimPlot(
        seurat_obj,
        reduction = "umap",
        label = F,
        repel = TRUE,
        pt.size = 0.5
      ) + NoLegend() +
        theme(aspect.ratio = 1, legend.position = "right")
      dev.off()




      meta_data <- seurat_obj@meta.data
      
      groundtruth_labels <- meta_data$Cancer_assignment
      
      table(groundtruth_labels)
      
      groundtruth_labels <- factor(
        ifelse(groundtruth_labels %in% c("Cancer"), 1, 0),
        levels = c(0, 1)
      )
      
      table(groundtruth_labels)
      
      seurat_obj$groundtruth_labels <- groundtruth_labels

      cairo_pdf(
        "/data1/liuchonghui/CTDP/figures of graduation/GSE182227_groundtruth_UMAP.pdf"
      )
      DimPlot(
        seurat_obj,
        group.by = "groundtruth_labels",
        reduction = "umap",
        label = F,
        repel = TRUE,
        pt.size = 0.5
      ) + NoLegend() +
        theme(aspect.ratio = 1, legend.position = "right")
      dev.off()


``````
