library(monocle3)
library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(dplyr)

curDate<-Sys.Date()

setwd("/home/flyhunter/Wang/output")

RNA.combined.norm<-readRDS('integ_OPC_ODC')

DefaultAssay(RNA.combined.norm)<-'RNA'

Idents(RNA.combined.norm)<-'integrated_snn_res.0.2'

targDir <- 'OPC_ODC/Monocle3/'

#Idents(RNA.combined.norm)<-'integrated_snn_res.0.3'

DimPlot(RNA.combined.norm, reduction = "umap", label = T, repel=T)

cds <- as.cell_data_set(RNA.combined.norm)


cds <- preprocess_cds(cds, num_dim = 85)

#cds <- preprocess_cds(cds)

#plot_pc_variance_explained(cds)


cds <- reduce_dimension(cds)

#plot_cells(cds)

plot_cells(cds, color_cells_by="integrated_snn_res.0.2")

colnames(cds@colData)

metadat<-data.frame(cds@colData)

cds <- align_cds(cds, num_dim = 85, alignment_group = "group")
cds <- reduce_dimension(cds)
conStr<-plot_cells(cds, color_cells_by="group", label_cell_groups=FALSE, cell_size=1)

ggsave(paste0(targDir, 'OPC_ODC_MonocClust_85PC_Root2_StressVsContr', curDate, '.jpeg'), plot = conStr, height = 12, width = 12, units = 'in', dpi = 300)

#cds <- cluster_cells(cds, resolution=1e-5)
cds <- cluster_cells(cds)
plot_cells(cds)

plot_cells(cds, color_cells_by = "partition")

q1<-cds@clusters$UMAP

clustID<-q1$clusters

cellNames<-names(clustID)
cellClusters<-as.character(clustID)
cellClustInf<-data.frame(cbind(cellNames, cellClusters))

### pseudotime
# learn trajectory
cds <- learn_graph(cds, use_partition = T)

regPlot<-monocle3::plot_cells(cds, group_label_size = 5, graph_label_size = 4)

regPlot

cds<-order_cells(cds, reduction_method = 'UMAP', root_cells = colnames(cds[, clusters(cds) =='2']))

psedPlot<-plot_cells(cds, group_label_size = 5, graph_label_size = 4, 
                     color_cells_by = 'pseudotime',
                     label_branch_points = F,
                     label_roots = F,
                     label_leaves = F, cell_size=1)+
  theme(text = element_text(size = 20))  

psedPlot

plotClust<-plot_cells(cds, group_label_size = 5, graph_label_size = 4, 
                      color_cells_by = 'cluster',
                      label_branch_points = F,
                      label_roots = F,
                      label_leaves = F, cell_size=1)+
  theme(text = element_text(size = 20))  

plotClust

ggsave(paste0(targDir, 'OPC_ODC_MonocClust_85PC_Root2_ColClust', curDate, '.jpeg'), plot = plotClust, height = 12, width = 12, units = 'in', dpi = 300)

cds$monocle3_pseudotime<- pseudotime(cds)

data.pseudo<- data.frame(colData(cds))

identical(row.names(data.pseudo), cellClustInf$cellNames)

data.pseudo$monocClust<-cellClustInf$cellClusters

psedBox<-ggplot(data.pseudo, aes(monocle3_pseudotime, cellClusters, fill= cellClusters))+
  geom_boxplot()+
  theme(text = element_text(size = 20)) 

psedBox

medTimes<-data.pseudo%>%
  group_by(monocClust)%>% 
  summarise(Mean=mean(monocle3_pseudotime), Median=median(monocle3_pseudotime))

ggsave(paste0(targDir,'OPC_ODC_MonocClust_85PC_Root2_', curDate, '.jpeg'), plot = psedPlot, height = 12, width = 12, units = 'in', dpi = 300)
ggsave(paste0(targDir, 'OPC_ODC_MonocClust_85PC_Root2_Box_', curDate, '.jpeg'), plot = psedBox, height = 12, width = 12, units = 'in', dpi = 300)
write.csv(medTimes, paste0(targDir, 'OPC_ODC_MonocClust_85PC_Root2_sum_', curDate, '.csv'), row.names = F)

monocUmap<-plot_cells(cds, label_branch_points = F,label_roots = F, cell_size=1, group_label_size = 4)+
  theme(text = element_text(size = 20)) 
ggsave(paste0(targDir, 'OPC_ODC_MonocClust_85PC_UMAP_', curDate, '.jpeg'), plot = monocUmap, height = 12, width = 12, units = 'in', dpi = 300)

oldLabel<-plot_cells(cds, color_cells_by="integrated_snn_res.0.2", label_roots = F, label_leaves = F, group_label_size = 4, cell_size=1)+
  theme(text = element_text(size = 20)) 
ggsave(paste0(targDir, 'OPC_ODC_MonocClust_85PC_UMAP_OldLabels', curDate, '.jpeg'), plot = oldLabel, height = 12, width = 12, units = 'in', dpi = 300)

# extract Umap Coordinates
traj.plot <- plot_cells(cds)
point.data <- ggplot_build(traj.plot)[["plot"]][["data"]]

umapDf<-point.data[, c('sample_name', 'data_dim_1', 'data_dim_2', 'CellName', 'group', 'cell_group')]

umapDf$scvelo_name <- NA

for (i in 1:nrow(umapDf)){
  if (umapDf$group[i]=='Control'){
    umapDf$scvelo_name[i]<-paste0('control:', gsub('-1', 'x', umapDf$CellName[i]))
  } else if (umapDf$group[i]=='Stress'){
    umapDf$scvelo_name[i]<-paste0('stress:', gsub('-1', 'x', umapDf$CellName[i]))
  }
}

scUmap<-umapDf[, c('scvelo_name', 'data_dim_1', 'data_dim_2')]
colnames(scUmap)[2:3]<-c('UMAP1', 'UMAP2')

write.csv(scUmap, 'OpcOdcInt_Mon3Clust85PC_UmapCoordScVel.csv', row.names = F)

write.csv(umapDf, 'OpcOdcInt_Mon3Clust85PC_scVelMetadata.csv', row.names = F)