library(Seurat)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(clustree)


DefaultAssay(RNA.combined.norm) <- "RNA"

DefaultAssay(RNA.combined.norm)
levels(RNA.combined.norm)

metadata<-data.frame(RNA.combined.norm@meta.data)

curDate<-Sys.Date()

genes<-c('Rbfox3', 'Snap25', 'Syt1', 'Fibcd1', 'Mpped1', 'Satb2', 'Cntn6', 'Kcnh5', 'Vwc2l', 'Nectin3',
         'Trps1', 'Glis3', 'Slc4a4', 'Gad1', 'Gad2', 'Plp1', 'Ptgds', 'Trf', 'Calcrl',
         'Cspg4', 'Vcan', 'Cx3cr1', 'Hexb', 'P2ry12', 'Ndnf', 'Reln')

clusters<-c('CA1', 'CA2', 'CA3', 'DG', 'GABA', 'ODC', 'OPC', 'MG','CR')

dotPlot<-DotPlot(object = RNA.combined.norm, features = genes, scale.max = 100, dot.scale = 16)+
  scale_colour_gradient2(low = "blue", mid = "grey", high = "red")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  theme(text = element_text(size = 24))+ # all text elements size
  theme(axis.text = element_text(size = 24))  # axes text size
  #scale_y_discrete(limits=rev)

ggsave(paste0('./dotPlotRenameClust0.3Res_', curDate, '.jpeg'),plot = dotPlot,  height = 12, width = 24, units = 'in', dpi = 300)

#RNA.combined.norm[["Annotations"]] <- Idents(object = RNA.combined.norm)

# identify markers per cluster

# all markers 
allMarkers <- FindAllMarkers(RNA.combined.norm , only.pos = F, min.pct = 0.1, logfc.threshold = 0.25, test.use = "MAST")

write.csv(allMarkers, paste0('allMarkersRenameClust0.3Res_',curDate, '.csv'), row.names = F)

allMarkersSign<-allMarkers[(allMarkers$p_val_adj<0.05),]

allMarkSignSum<-data.frame(table(allMarkersSign$cluster))
colnames(allMarkSignSum)<-c('Cluster', 'Significant_markers_number')
write.csv(allMarkSignSum, 'sign_all_markers_sum_RenameClust0.3Res.csv', row.names = F)

topMarkers<-character()

extractTopMarkers<-function(x){
  topMarkersDf<-data.frame(matrix(nrow = 0, ncol = 0))
  for ( i in clusters){
    df<-x[(x$cluster==i),]
    df$AbsLog2FC<-abs(df$avg_log2FC)
    dfSort<-df[order(df$AbsLog2FC, decreasing = T),]
    top10<-head(dfSort, 10)
    topMarkersDf<-rbind(topMarkersDf, top10)
  }
  return(topMarkersDf)
}
top10MarkersDf<-extractTopMarkersPval(allMarkersSign)

topMarkers<-unique(top10MarkersDf$gene)

allHeat<-DoHeatmap(RNA.combined.norm, features = topMarkers) + NoLegend()

ggsave(paste0('./heatMap_top10Genes_RenameClust0.3Res_', curDate, '.jpeg'), plot = allHeat, height = 16, width = 18, units = 'in', dpi = 300)


allHeat<-DoHeatmap(RNA.combined.norm, features = topMarkers)

ggsave(paste0('./heatMapLeg_top10Genes_RenameClust0.3Res_', curDate, '.jpeg'), plot = allHeat, height = 20, width = 20, units = 'in', dpi = 300)

rm(allHeat, allMarkers, allMarkersSign, top10MarkersDf, topMarkers)
# only positive markers
posMarkers <- FindAllMarkers(RNA.combined.norm , only.pos = T, min.pct = 0.1, logfc.threshold = 0.25, test.use = "MAST")

write.csv(posMarkers, paste0('Positive_MarkersRenameClust0.3Res_',curDate, '.csv'), row.names = F)

posMarkersSign<-posMarkers[(posMarkers$p_val_adj<0.05),]

top10MarkersDf<-extractTopMarkers(posMarkersSign)

topMarkers<-unique(top10MarkersDf$gene)

posHeat<-DoHeatmap(RNA.combined.norm, features = topMarkers) + NoLegend()

ggsave(paste0('./heatMap_top10PositiveGenes_RenameClust0.3Res_', curDate, '.jpeg'), plot = posHeat, height = 16, width = 18, units = 'in', dpi = 300)


posHeat<-DoHeatmap(RNA.combined.norm, features = topMarkers)

ggsave(paste0('./heatMapLeg_top10PositiveGenes_RenameClust0.3Res_', curDate, '.jpeg'), plot = posHeat, height = 20, width = 20, units = 'in', dpi = 300)

rm(posHeat, posMarkers, posMarkersSign, top10MarkersDf, topMarkers)

###
