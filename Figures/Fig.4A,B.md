## Fig4.A,B
### Download and preprocess the single-cell data of melanoma
```r
  library(Seurat)
  library(dplyr)
  expr_file <- "/data/chliu/Supervised_Phen_subpopulations_project/data/GSE120575/GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt"
  meta_file <- "/data/chliu/Supervised_Phen_subpopulations_project/data/GSE120575/GSE120575_patient_ID_single_cells.txt"
  expr_data <- read.delim(expr_file, header = TRUE, sep = "\t", row.names = 1)
  expr_data[1:5,1:5]
  expr_data <- expr_data[-1, -1]
  expr_data[1:5,1:5]
  dim(expr_data)
  expr_data <- as.data.frame(expr_data)
  meta_data <- read.table(meta_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, fill = TRUE) 
  dim(meta_data)
  meta_data[1:5,1:5]
  meta_data <-meta_data[,which(colnames(meta_data) %in% c("title","characteristics..response"))]
  meta_data[1:5,]
  colnames(meta_data) <- c("cell_id", "response")
  rownames(meta_data) <- meta_data$cell_id
  meta_data[1:5,]
  all(colnames(expr_data) == rownames(meta_data))
  rownames(meta_data)=colnames(expr_data)
  meta_data$cell_id =colnames(expr_data)
  all(colnames(expr_data) == rownames(meta_data))
  GSE120575_sc_data <- CreateSeuratObject(counts = expr_data, meta.data = meta_data, min.cells = 3, min.features = 200)
  print(GSE120575_sc_data)
  GSE120575_sc_data[["percent.mt"]] <- PercentageFeatureSet(GSE120575_sc_data, pattern = "^MT-")
  GSE120575_sc_data <- subset(GSE120575_sc_data, subset = percent.mt < 20)
  dim(GSE120575_sc_data@assays[["RNA"]]@counts)
  GSE120575_sc_data <- NormalizeData(GSE120575_sc_data, normalization.method = "LogNormalize", scale.factor = 10000)
  GSE120575_sc_data <- ScaleData(GSE120575_sc_data)
  GSE120575_sc_data <- FindVariableFeatures(GSE120575_sc_data, selection.method = "vst", nfeatures = 2000)
  GSE120575_sc_data <- RunPCA(GSE120575_sc_data, features = VariableFeatures(object = GSE120575_sc_data))
  GSE120575_sc_data <- FindNeighbors(GSE120575_sc_data, dims = 1:10)
  GSE120575_sc_data <- FindClusters(GSE120575_sc_data)
  GSE120575_sc_data <- RunUMAP(GSE120575_sc_data, dims = 1:10)
  save(GSE120575_sc_data,file="/data/chliu/Supervised_Phen_subpopulations_project/test/GSE120575_sc_data.rdata")
  
``````
### Downsampling of single-cell data in melanoma cases
```r
  load("/data/chliu/Supervised_Phen_subpopulations_project/test/GSE120575_sc_data.rdata")
  table(GSE120575_sc_data@meta.data$response)
  library(dplyr)
  library(Seurat)
  meta_data <- GSE120575_sc_data@meta.data
  target_count <- sum(meta_data$response == "Responder")
  non_responder_counts <- meta_data %>%
    filter(response == "Non-responder") %>%
    group_by(seurat_clusters) %>%
    summarise(cluster_count = n())
  non_responder_counts <- non_responder_counts %>%
    mutate(target_sample_size = round(cluster_count / sum(cluster_count) * target_count))
  set.seed(123) 
  sampled_non_responders <- data.frame()
  for (i in 1:nrow(non_responder_counts)) {
    cluster_id <- non_responder_counts$seurat_clusters[i]
    sample_size <- non_responder_counts$target_sample_size[i]
    
    cells_in_cluster <- meta_data %>%
      filter(response == "Non-responder", seurat_clusters == cluster_id)
    
    sampled_cells <- if (nrow(cells_in_cluster) <= sample_size) {
      cells_in_cluster
    } else {
      cells_in_cluster[sample(1:nrow(cells_in_cluster), sample_size), ]
    }
    
    sampled_non_responders <- rbind(sampled_non_responders, sampled_cells)
  }
  responders <- meta_data %>%
    filter(response == "Responder")
  downsampled_meta_data <- bind_rows(sampled_non_responders, responders)
  GSE120575_sc_data_downsampled <- subset(GSE120575_sc_dalta, cells = downsampled_meta_data$cell_id)
  table(downsampled_meta_data$response)
  dim(GSE120575_sc_data_downsampled)
  save(GSE120575_sc_data_downsampled,file = "/data/chliu/Supervised_Phen_subpopulations_project/test/GSE120575_sc_data_downsampled.rdata")
   
  
``````
### UAMP plots of single-cell data for cases of melanoma
```r
  load("/data/chliu/Supervised_Phen_subpopulations_project/test/GSE120575_sc_data_downsampled.rdata")
  library(Seurat)
  library(ggpubr)
  library(ggplot2)
  cairo_pdf(paste0("/data/chliu/Supervised_Phen_subpopulations_project/figures of graduation/GSE120575_sc_data_downsampled_UMAP.pdf"))
  DimPlot(GSE120575_sc_data_downsampled,  reduction = "umap", label = T, repel = TRUE,pt.size = 0.5) + NoLegend()+theme(aspect.ratio=1, legend.position = "right")
  dev.off()
  
  cairo_pdf(paste0("/data/chliu/Supervised_Phen_subpopulations_project/figures of graduation/GSE120575_sc_data_downsampled_source_phen_UMAP.pdf"))
  DimPlot(GSE120575_sc_data_downsampled, group.by = "response", reduction = "umap", label = F, repel = TRUE,pt.size = 0.5) + NoLegend()+theme(aspect.ratio=1, legend.position = "right")
  dev.off()  
  
``````