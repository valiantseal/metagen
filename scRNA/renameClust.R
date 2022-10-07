library(Seurat)
library(dplyr)
library(ggplot2)




setwd("C:/Users/abomb/Projects/AU/wang/SingleCellStress/output")

RNA.combined.norm<-readRDS('intergratedRNA')


DefaultAssay(RNA.combined.norm) <- "integrated"
DefaultAssay(RNA.combined.norm)
Idents(RNA.combined.norm)<-'integrated_snn_res.0.3'

metadata<-data.frame(RNA.combined.norm@meta.data)

currentClustNames<-levels(RNA.combined.norm)
newClustNames<-currentClustNames

for ( i in 1:length(newClustNames)){
  if (newClustNames[i] == 0 | newClustNames[i] == 1 | newClustNames[i] == 3){
    newClustNames[i]<-'DG'
  } else if ( newClustNames[i] == 2){
    newClustNames[i]<-'CA1'
  } else if ( newClustNames[i] == 4 | newClustNames[i] == 10 ) {
    newClustNames[i]<-'ODC'
  } else if ( newClustNames[i] == 5 ) {
    newClustNames[i]<-'OPC'
  } else if ( newClustNames[i] == 6 | newClustNames[i] == 12 ) {
    newClustNames[i]<-'GABA'
  } else if ( newClustNames[i] == 8 ) {
    newClustNames[i]<-'CA3'
  } else if ( newClustNames[i] == 7 ) {
    newClustNames[i]<-'MG'
  } else if ( newClustNames[i] == 9 | newClustNames[i] == 11 | newClustNames[i] == 14 ) {
    newClustNames[i]<-'CA2'
  } else if ( newClustNames[i] == 13 ) {
    newClustNames[i]<-'CR'
  }
}

names(newClustNames) <- levels(RNA.combined.norm)
RNA.combined.norm <- RenameIdents(RNA.combined.norm, newClustNames)

levels(x =RNA.combined.norm) <- c('CA1', 'CA2', 'CA3', 'DG', 'GABA', 'ODC', 'OPC', 'MG','CR')

curDate<-Sys.Date()


# make plots
DimPlot(RNA.combined.norm, reduction = "umap", label = F)
ggsave( paste0('./combCondUmapRenameClust0.3Res_', curDate, '.jpeg'), height = 6, width = 10, units = 'in', dpi = 300)

control<-subset(x = RNA.combined.norm, subset = group == "Control")
DimPlot(control, reduction = "umap", label = F)
ggsave( paste0('./controlUmapRenameClust0.3Res_', curDate, '.jpeg'), height = 6, width = 10, units = 'in', dpi = 300)

stress<-subset(x = RNA.combined.norm, subset = group == "Stress")
DimPlot(stress, reduction = "umap", label = F)
ggsave( paste0('./stressUmapRenameClust0.3Res_', curDate, '.jpeg'), height = 6, width = 10, units = 'in', dpi = 300)


# count cells per cluster
nrow(control@meta.data)
nrow(stress@meta.data)

RNA.combined.norm[["Annotations"]] <- Idents(object = RNA.combined.norm)
metadata<-data.frame(RNA.combined.norm@meta.data)


cellsPerClust<-data.frame(table(metadata$Annotations, metadata$group))
colnames(cellsPerClust)<-c('Cluster', 'Group', 'Cells_number')
write.csv(cellsPerClust, paste0('./Cells_per_cluster_group_res0.3_', curDate, '.csv'), row.names = F)


rm(control, stress, metadata)

