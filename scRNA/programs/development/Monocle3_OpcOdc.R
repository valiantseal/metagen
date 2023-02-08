library(monocle3)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(dplyr)

curDate<-Sys.Date()

setwd("C:/Users/abomb/Projects/AU/wang/SingleCellStress/output")

RNA.combined.norm<-readRDS('integ_OPC_ODC')

DefaultAssay(RNA.combined.norm)<-'RNA'

Idents(RNA.combined.norm)<-'integrated_snn_res.0.2'

#Idents(RNA.combined.norm)<-'integrated_snn_res.0.3'

DimPlot(RNA.combined.norm, reduction = "umap", label = T, repel=T)

cds <- as.cell_data_set(RNA.combined.norm)

colData(cds)

fData(cds)$gene_short_name <- rownames(fData(cds))

recreate.partition<- c(rep(1, length(cds@colData@rownames)))

names(recreate.partition) <- cds@colData@rownames

recreate.partition <- as.factor(recreate.partition)

cds@clusters$UMAP$partitions <- recreate.partition

list_cluster<- RNA.combined.norm@active.ident

cds@clusters$UMAP$clusters<-list_cluster

cds@int_colData@listData$reducedDims$UMAP<-RNA.combined.norm$umap@cell.embeddings

monocle3::plot_cells(cds)

#metaSel$New_ID<-rownames(metaSel)

metaOO<-RNA.combined.norm@meta.data
metaOO$New_ID<-rownames(metaOO)

#write.csv(metaOO, 'OPC_ODC_combMeta.csv', row.names = F)

#combMeta<-plyr::join(metaOO, metaSel, by='New_ID', type='left', match='all')

# learn trajectory
cds <- learn_graph(cds)

regPlot<-monocle3::plot_cells(cds, group_label_size = 5, graph_label_size = 4)

#ggsave(paste0('OPC_ODC_combinedMonoc_NotOrd_', curDate, '.jpeg'), plot = regPlot, height = 12, width = 12, units = 'in', dpi = 300)

cds<-order_cells(cds, reduction_method = 'UMAP', root_cells = colnames(cds[, clusters(cds) =='1']))

psedPlot<-plot_cells(cds, group_label_size = 5, graph_label_size = 4, 
                     color_cells_by = 'pseudotime',
                     label_branch_points = F,
                     label_roots = F,
                     label_leaves = F)

psedPlot
#ggsave(paste0('OPC_ODC_Res0.3_Root1_', curDate, '.jpeg'), plot = psedPlot, height = 12, width = 12, units = 'in', dpi = 300)

cds$monocle3_pseudotime<- pseudotime(cds)

data.pseudo<- data.frame(colData(cds))

psedBox<-ggplot(data.pseudo, aes(monocle3_pseudotime, integrated_snn_res.0.2, fill= integrated_snn_res.0.2))+
geom_boxplot()

psedBox

ggsave(paste0('OPC_ODC_Res0.3_Root1_Box_', curDate, '.jpeg'), plot = psedBox, height = 12, width = 12, units = 'in', dpi = 300)

medTimes<-data.pseudo%>%
  group_by(integrated_snn_res.0.2)%>% 
  summarise(Mean=mean(monocle3_pseudotime), Median=median(monocle3_pseudotime))

#write.csv(medTimes, paste0('OPC_ODC_Res0.2_Root1_sum_', curDate, '.csv'), row.names = F)

