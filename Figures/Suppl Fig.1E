### The PENCIL algorithm is executed
```r
        library(Seurat)
        library(SeuratDisk)
        load("/data1/liuchonghui/CTDP/test/GSE182227_sc_data.rdata")
        
        meta_data <- GSE182227_sc_data@meta.data
        cell.labels <- ifelse(
          meta_data$Sequenced_sample %in% c(
            "OP33",
            "OP34",
            "OP35BOT"
          ),
          "1",
          "0"
        )
        table(meta_data$Sequenced_sample)
        table(cell.labels)
        GSE182227_sc_data <- AddMetaData(GSE182227_sc_data, metadata = cell.labels,col.name = "source_phen")
        GSE182227_sc_data@meta.data[1:5,10:15]
        

        library(Seurat)
        library(SingleCellExperiment)
        library(zellkonverter)
  
        scale_mat <- GetAssayData(GSE182227_sc_data, assay = "RNA", layer = "scale.data")
        dim(scale_mat) 
        
        sce <- as.SingleCellExperiment(GSE182227_sc_data, assay = "RNA")

        assay(sce, "scaled") <- scale_mat

        writeH5AD(
          sce,
          file = "/data1/liuchonghui/CTDP/test/GSE182227_sc_data_scaled.h5ad",
          X_name = "scaled"
        )
``````
```python
import os
        os.environ['CUDA_VISIBLE_DEVICES'] = '-1'
        
        from pencil import *
        import scanpy as sc
        from matplotlib.colors import LinearSegmentedColormap
        import warnings
        warnings.filterwarnings('ignore') #to ignore warnings on the output.
        new_directory = "/data1/liuchonghui/CTDP/test/"
        os.chdir(new_directory)
        print("当前工作目录已修改为:", os.getcwd())
        
        
        adata = sc.read_h5ad('/data1/liuchonghui/CTDP/test/GSE182227_sc_data_scaled.h5ad')
        print(adata)
        
        num_classes = len(set(adata.obs.source_phen.values))
        labels_raw = adata.obs.source_phen
        labels_raw = pd.Categorical(labels_raw)
        data, labels = adata.X.copy(), labels_raw.codes
        class_names = list(labels_raw.categories)
        
        print(data.shape)
        print(labels.shape)
        print(class_names)
        
        
        data_name = 'GSE182227_sc_data_scaled'
        expr_id = '0.0.1'
        mode = 'multi-classification'
        pencil = Pencil(mode, select_genes=True, seed=1234, data_name=data_name, expr_id=expr_id, mlflow_record=True)
        with mlflow.start_run():
          pred, confidence = pencil.fit_transform(
            data, labels, 
            class_names=class_names,
            plot_show=True
          )
        w = pencil.gene_weights(plot=True)
        plt.close()
        
        
        pred_labels = np.array(labels_raw.categories[pred])
        pred_labels=pred_labels.astype(str)
        pred_labels[confidence < 0] = 'Rejected'
        np.savetxt("GSE182227_sc_data_scaled_PENCILpred_labels.csv", pred_labels, delimiter=",", fmt="%s")
``````
### Generate an UMAP plot
```r
      load("/data1/liuchonghui/CTDP/test/GSE182227_sc_data.rdata")
      seurat_obj <- GSE182227_sc_data

      predicted_labels <- read.csv(
        "/data1/liuchonghui/CTDP/test/GSE182227_sc_data_scaled_PENCILpred_labels.csv",
        header = FALSE
      )
      
      table(predicted_labels)
      
      pred_labels <- as.data.frame(
        ifelse(predicted_labels == "1", 1, 0)
      )
      table(pred_labels)
      
      colnames(pred_labels) <- "PENCIL_pred_labels"
      rownames(pred_labels) <- rownames(seurat_obj@meta.data)
      
      seurat_obj <- AddMetaData(seurat_obj, metadata = pred_labels)
      colnames(seurat_obj@meta.data)

      library(Seurat)
      library(ggpubr)
      library(ggplot2)
      cairo_pdf(
        "/data1/liuchonghui/CTDP/figures of graduation/GSE182227_PENCIL_pred_UMAP.pdf"
      )
      DimPlot(
        seurat_obj,
        group.by = "PENCIL_pred_labels",
        reduction = "umap",
        label = F,
        repel = TRUE,
        pt.size = 0.5
      ) + NoLegend() +
        theme(aspect.ratio = 1, legend.position = "right")
      dev.off()
      
``````
