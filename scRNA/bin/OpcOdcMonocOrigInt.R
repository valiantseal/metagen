library(monocle3)
library(Seurat)
library(SeuratWrappers)
library(gridExtra)

curDate<-Sys.Date()

setwd("C:/Users/abomb/Projects/AU/wang/SingleCellStress/output")

source('../../programs/renameClusters.R')

selCells<-subset(x = RNA.combined.norm, idents = c("OPC", "ODC"))

DefaultAssay(selCells)<-'integrated'

rm(RNA.combined.norm)

RNA.combined.norm<-selCells

RNA.combined.norm <- RunPCA(RNA.combined.norm, verbose = FALSE)

###
#ElbowPlot(RNA.combined.norm, ndims = 50)
#ggsave("ElbowPlot_DupFiltr.jpeg", width = 16, height = 10)

# function borrowed from Sergei Bombin
fgetNPCs <- function(obj,MIN_CUM_SD){
  sds <- Stdev(obj,reduction="pca")
  cumsds <- 0
  for(i in 1:length(sds)){
    cumsds <- cumsds + sds[i]
    if(cumsds >= MIN_CUM_SD){
      return(i)
    }
  }
}

#minPCA<-fgetNPCs(obj=RNA.combined.norm, MIN_CUM_SD=95)

#RNA.combined.norm <- JackStraw(RNA.combined.norm, num.replicate = 100, dims = 30)
#RNA.combined.norm <- ScoreJackStraw(RNA.combined.norm, dims = 1:30)
#JackStrawPlot(RNA.combined.norm, dims = 1:30)

RNA.combined.norm <- FindNeighbors(RNA.combined.norm, dims = 1:15)
RNA.combined.norm <- RunUMAP(RNA.combined.norm, reduction = "pca", dims = 1:15)
resols<-c(0,0.05,0.08,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)
RNA.combined.norm <- FindClusters(RNA.combined.norm, resolution = resols)

i<-0.2
current_resol<-paste0("integrated_snn_res.", i)
Idents(RNA.combined.norm)<-current_resol
p1 <- DimPlot(RNA.combined.norm, reduction = "umap", group.by = "group", label = F)+ggtitle(current_resol)
p2 <- DimPlot(RNA.combined.norm, reduction = "umap", label = T, repel=T)+ggtitle(current_resol)
p3<-grid.arrange(p1,p2, ncol=2, nrow=1)


ggsave(paste0('OPC_ODC_umapClustering_OrigInt_',current_resol, "_", curDate, ".jpeg"), plot=p3, height = 6, width = 10, units = 'in', dpi = 300)
# monocle 3
DefaultAssay(RNA.combined.norm)<-'RNA'

Idents(RNA.combined.norm)<-'integrated_snn_res.0.2'

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

cds<-order_cells(cds, reduction_method = 'UMAP', root_cells = colnames(cds[, clusters(cds) =='3']))

psedPlot<-plot_cells(cds, group_label_size = 5, graph_label_size = 4, 
                     color_cells_by = 'pseudotime',
                     label_branch_points = F,
                     label_roots = F,
                     label_leaves = F)

psedPlot
ggsave(paste0('OPC_ODC_Res0.3_Root3_', curDate, '.jpeg'), plot = psedPlot, height = 12, width = 12, units = 'in', dpi = 300)

cds$monocle3_pseudotime<- pseudotime(cds)

data.pseudo<- data.frame(colData(cds))

psedBox<-ggplot(data.pseudo, aes(monocle3_pseudotime, integrated_snn_res.0.2, fill= integrated_snn_res.0.2))+
  geom_boxplot()

psedBox

ggsave(paste0('OPC_ODC_Res0.2_OrigInt_Root3_Box_', curDate, '.jpeg'), plot = psedBox, height = 12, width = 12, units = 'in', dpi = 300)

medTimes<-data.pseudo%>%
  group_by(integrated_snn_res.0.2)%>% 
  summarise(Mean=mean(monocle3_pseudotime), Median=median(monocle3_pseudotime))

write.csv(medTimes, paste0('OPC_ODC_Res0.2_OrigInt_Root3_sum_', curDate, '.csv'), row.names = F)
