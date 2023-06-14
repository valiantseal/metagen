library(clusterProfiler)
library(enrichplot)
library(org.Mm.eg.db)

targDir = './Paper_figs/Fig5/'

markers = read.csv(paste0(targDir, "DiffExprRNAMonoc3ClustAdjP_2023-06-08.csv"))
clusters = unique(markers$Cluster)
curDate = Sys.Date()

dfFilt<-markers[(markers$fdr_p < 0.05),]
entrez <- as.data.frame(mapIds(org.Mm.eg.db, keys = dfFilt$Genes, keytype = "SYMBOL", column="ENTREZID"))
colnames(entrez)<-'entrez'
dfEntr<-cbind(dfFilt, entrez)

allEnrichedKegg<-function(dataTab, clusters, targDir, curDate){
  combDat = data.frame(matrix(nrow = 0, ncol = 0))
  for ( cluster in clusters){
    dfSel<-dataTab[(dataTab$Cluster==cluster),]
    if (nrow(dfSel) > 0) {
      genes<-unique(dfSel$entrez)
      
      ego2 <- enrichKEGG(gene = genes, organism = 'mmu',
                         pvalueCutoff = 0.05)
      
      ego2Res<-ego2@result
      ego2Res$Cluster = cluster
      
      dir.create(targDir, showWarnings = FALSE, recursive = T)
      
      tablePath<-paste0(targDir, cluster, '_', curDate, '.csv')
      write.csv(ego2Res, tablePath, row.names = F)
      combDat = rbind(combDat, ego2Res)
      
      if (nrow(ego2Res[(ego2Res$p.adjust < 0.05),]) > 0){
        plot1<-barplot(ego2, showCategory=20)
        plotPath<-paste0(targDir, cluster, '_20Terms_', curDate, '.jpeg')
        
        png(filename = plotPath, height = 16, width = 20, units = 'in', res = 300)
        print(plot1)
        dev.off()
        
      }
    }
    
  }
  return(combDat)
}

keggAll = allEnrichedKegg(dataTab=dfEntr, clusters = clusters, targDir = './Paper_figs/Fig5/ClustProf/Kegg_all/', curDate = curDate)
write.csv(keggAll, paste0(targDir, "AllClust_AllGenes_2023-06-12.csv"), row.names = F)

# only positive 
dfEntr<-cbind(dfFilt, entrez)
dfEntr = dfEntr[(dfEntr$avg_log2FC > 0),]
keggPos = allEnrichedKegg(dataTab=dfEntr, clusters = clusters, targDir = './Paper_figs/Fig5/ClustProf/Kegg_positive/', curDate = curDate)
write.csv(keggPos, paste0(targDir, "AllClust_PosGenes_2023-06-12.csv"), row.names = F)

# only negative 
dfEntr<-cbind(dfFilt, entrez)
dfEntr = dfEntr[(dfEntr$avg_log2FC < 0),]
keggNeg = allEnrichedKegg(dataTab=dfEntr, clusters = clusters, targDir = './Paper_figs/Fig5/ClustProf/Kegg_negative/', curDate = curDate)
write.csv(keggNeg, paste0(targDir, "AllClust_NegGenes_2023-06-12.csv"), row.names = F)
