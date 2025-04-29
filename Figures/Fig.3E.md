## Fig3.E
### The PENCIL algorithm is executed
```r
  library(Seurat)
  library(SeuratDisk)
  library(caret)
  load("/data/chliu/Supervised_Phen_subpopulations_project/test/PBMC_Travaglini_SC_sample3clusters_list.rdata")
  process_single_data <- function(pbmc, index) {
    SaveH5Seurat(pbmc, filename = paste0("/data/chliu/Supervised_Phen_subpopulations_project/test/PBMC_Travaglini_SC_3clusterslabeled_", index,".h5Seurat")) 
    Convert(paste0("/data/chliu/Supervised_Phen_subpopulations_project/test/PBMC_Travaglini_SC_3clusterslabeled_", index,".h5Seurat"), dest = "h5ad")
  }
  lapply(seq_along(PBMC_Travaglini_SC_sample3clusters_list), function(i) {
    process_single_data(PBMC_Travaglini_SC_sample3clusters_list[[i]], i-1)
  })
``````
```python
from pencil import *
import scanpy as sc
import numpy as np
import pandas as pd
from sklearn.metrics import f1_score
import os
import mlflow
import warnings
warnings.filterwarnings('ignore')
os.environ['CUDA_VISIBLE_DEVICES'] = '3'

# 设置工作目录
new_directory = "/home/liuchonghui/Supervised_Phen_subpopulations_project/test/"
os.chdir(new_directory)
print("当前工作目录已修改为:", os.getcwd())


# 定义处理单个对象的函数
def process_single_h5ad(index):
    # 读取.h5ad文件
    adata = sc.read_h5ad(f'PBMC_Travaglini_SC_3clusterslabeled_{index}.h5ad')
    
    # 准备数据
    labels_raw = adata.obs['source_phen']
    labels_raw = pd.Categorical(labels_raw)
    data, labels = adata.X.copy(), labels_raw.codes
    class_names = list(labels_raw.categories)
    
    # 执行PENCIL算法
    data_name = f'PBMC_Travaglini_SC_3clusterslabeled_{index}'
    expr_id = '0.0.1'
    mode = 'multi-classification'
    pencil = Pencil(mode, select_genes=True, seed=1234, data_name=data_name, expr_id=expr_id, mlflow_record=True)
    
    with mlflow.start_run():
        pred, confidence = pencil.fit_transform(data, labels, class_names=class_names, plot_show=True)
        pred_labels = np.array(labels_raw.categories[pred]).astype(str)
        pred_labels[confidence < 0] = 'Rejected'
        
    # 保存预测标签
    np.savetxt(f"PBMC_Travaglini_SC_3clusterslabeled_{index}_PENCILpred_labels.csv", pred_labels, delimiter=",", fmt="%s")
    
    # 返回预测标签
    return pred_labels

# 对列表中的每个数据对象执行处理
for i in range(10):
    pred_labels = process_single_h5ad(i)


``````
### Generate an UMAP plot
```r
  load("/data/chliu/Supervised_Phen_subpopulations_project/test/PBMC_Travaglini_SC_sample3clusters_list.rdata")
  i = 5
  print(names(PBMC_Travaglini_SC_sample3clusters_list)[i])
  seurat_obj <- PBMC_Travaglini_SC_sample3clusters_list[[i]]
  
  predicted_labels <- read.csv(paste0("/data/chliu/Supervised_Phen_subpopulations_project/test/PBMC_Travaglini_SC_3clusterslabeled_", i-1, "_PENCILpred_labels.csv"), header = FALSE)
  table(predicted_labels)
  pred_labels <- as.data.frame(ifelse(predicted_labels == "1.0", 1, 0))
  table(pred_labels)
  colnames(pred_labels) <- "PENCIL_pred_labels"
  rownames(pred_labels) <- rownames(seurat_obj@meta.data)
  seurat_obj <- AddMetaData(seurat_obj, metadata = pred_labels)
  seurat_obj@meta.data[1:5,16:ncol(seurat_obj@meta.data)]
  
  cairo_pdf(paste0("/data/chliu/Supervised_Phen_subpopulations_project/figures of graduation/",names(PBMC_Travaglini_SC_sample3clusters_list)[i],"_PENCIL_pred_UMAP.pdf"))
  DimPlot(seurat_obj, group.by = "PENCIL_pred_labels", reduction = "umap", label = F, repel = TRUE,pt.size = 0.5) + NoLegend()+theme(aspect.ratio=1, legend.position = "right")
  dev.off()
  
``````