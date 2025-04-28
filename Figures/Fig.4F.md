## Fig4.F
```r
 library(SCENIC)
  setwd("/data/chliu/Supervised_Phen_subpopulations_project/test/immutherapy_response_SCENIC.result/")
  
  load("/data/chliu/Supervised_Phen_subpopulations_project/test/GSE120575_sc_data_downsampled.rdata")
  exprMat <- as.matrix(Seurat::GetAssayData(object = GSE120575_sc_data_downsampled, slot = "data"));gc()
  exprMat[1:5,1:5]
  load("/data/chliu/Supervised_Phen_subpopulations_project/test/GSE120575_logistic_result.rdata")
  logistic_result <- GSE120575_logistic_result
  significant_clusters <- rownames(logistic_result[logistic_result$FDR < 0.05, ])
  print(significant_clusters)
  cellInfo <- as.data.frame(
    ifelse(GSE120575_sc_data_downsampled@meta.data$seurat_clusters %in% as.numeric(gsub("cluster", "", significant_clusters)), 1, 0)
  )
  table(cellInfo)
  colnames(cellInfo) <- "CellType"
  table(cellInfo[,1])
  
  scenicOptions <- initializeScenic(org="hgnc", dbDir="/data/chliu/PACSI_project/input_data", nCores=100)
  saveRDS(scenicOptions, file="int/scenicOptions.Rds")
  
  ### Co-expression network
  genesKept <- geneFiltering(exprMat, scenicOptions)
  exprMat_filtered <- exprMat[genesKept, ]
  exprMat_filtered[1:4,1:4]
  runCorrelation(exprMat_filtered, scenicOptions)
  exprMat_filtered_log <- exprMat_filtered
  runGenie3(exprMat_filtered_log, scenicOptions)
  
  exprMat_log <- exprMat
  scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"]
  scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
  library(data.table)
  source("/data/chliu/PACSI_project/input_data/runSCENIC_2_createRegulons.R")
  scenicOptions <- runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top5perTarget")) # Toy run settings
  library(foreach)
  scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log)
  
  
  
  regulonAUC <- readRDS(file="/data/chliu/Supervised_Phen_subpopulations_project/test/immutherapy_response_SCENIC.result/int/3.4_regulonAUC.Rds")
  library(AUCell)
  library(SCENIC)
  load("/data/chliu/Supervised_Phen_subpopulations_project/test/GSE120575_logistic_result.rdata")
  load("/data/chliu/Supervised_Phen_subpopulations_project/test/GSE120575_sc_data_downsampled.rdata")
  logistic_result <- GSE120575_logistic_result
  significant_clusters <- rownames(logistic_result[logistic_result$FDR < 0.05, ])
  print(significant_clusters)
  cellInfo <- as.data.frame(
    ifelse(GSE120575_sc_data_downsampled@meta.data$seurat_clusters %in% as.numeric(gsub("cluster", "", significant_clusters)), 1, 0)
  )
  colnames(cellInfo) <- "CellType"
  rss <- calcRSS(AUC=getAUC(regulonAUC), cellAnnotation=as.character(cellInfo$CellType))
  rssPlot <- plotRSS(rss)
  plotly::ggplotly(rssPlot$plot)
  options(repr.plot.width=5, repr.plot.height=5) # To set the figure size in Jupyter
  cairo_pdf(file="/data/chliu/Supervised_Phen_subpopulations_project/figures of graduation/immutherapy_response_SCENIC.pdf")
  plotRSS_oneSet(rss, setName = 1)
  dev.off() 
  
  
``````
