## Fig6.E
```r
  load("D://Supervised_Phen_subpopulations_project/test/Livercirrhosis_Cirrhotic_associatedcells_sig.rdata")
  data <- subset(Cirrhotic_associatedcells_sig, p_val_adj < 0.05 & avg_log2FC > 0.5)
  #symbol.list <- rownames(data)
  library(clusterProfiler);library(org.Hs.eg.db)
  gene.df <- bitr(rownames(data),fromType = "SYMBOL",
                  toType = c("ENTREZID"),OrgDb = org.Hs.eg.db)
  Reactome_enrich.result <- ReactomePA::enrichPathway(gene=gene.df$ENTREZID, pvalueCutoff = 0.05, readable=TRUE)
  save(Reactome_enrich.result, file = "D://Supervised_Phen_subpopulations_project/test/Livercirrhosis_Cirrhotic_Reactome_enrich.result.rdata")
  
  load("D://Supervised_Phen_subpopulations_project/test/Livercirrhosis_Cirrhotic_Reactome_enrich.result.rdata")
  library(ggplot2)
  cairo_pdf("D:/Supervised_Phen_subpopulations_project/figures of graduation/Livercirrhosis_Cirrhotic_reactome.barplot.pdf",width = 10)
  ggplot2::ggplot(Reactome_enrich.result, aes(
    y = -log10(Reactome_enrich.result@result$p.adjust[1:5]),
    x = reorder(Reactome_enrich.result@result$Description[1:5], -log10(Reactome_enrich.result@result$p.adjust[1:5]))
  )) + 
    geom_bar(stat = "identity", position = "dodge") + 
    coord_flip() + 
    theme_bw()
  dev.off()
``````
