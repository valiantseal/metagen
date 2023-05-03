library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(MAST)

curDate<-Sys.Date()

setwd("/home/flyhunter/Wang/output")

targDir = './Paper_figs/Fig3/'

dir.create(targDir, recursive = T, showWarnings = F)

RNA.combined.norm<-readRDS('integ_OPC_ODC')

DimPlot(RNA.combined.norm, group.by = "MonocClust")

RNA.combined.norm$newMonocClust = RNA.combined.norm$MonocClust

RNA.combined.norm$newMonocClust[RNA.combined.norm$newMonocClust == 4] = 3

RNA.combined.norm$newMonocClust[RNA.combined.norm$newMonocClust == 1] = "ODC"
RNA.combined.norm$newMonocClust[RNA.combined.norm$newMonocClust == 2] = "OPC"
RNA.combined.norm$newMonocClust[RNA.combined.norm$newMonocClust == 3] = "Intermideate"

Idents(RNA.combined.norm) = RNA.combined.norm$newMonocClust

levels(x =RNA.combined.norm) <- c('OPC', 'Intermideate', 'ODC')

DimPlot(RNA.combined.norm, group.by = "newMonocClust")

groupMarkers = read.csv("OPC_ODC/StressVsContr/RNA_DE/DiffExprMonoc3Clust_AllClust_ContrVsStress_2023-04-10.csv")

# find top genes
getTopMarkers = function(df, topNumb) {
  clusters = unique(df$Cell_Type)
  markers = character()
  dfPct = df[(df$pct.1 >0.25) | (df$pct.2 >0.25), ]
  for (cluster in clusters) {
    dfSub =  dfPct[(dfPct[['Cell_Type']] == cluster),]
    dfSub$AbsLog = abs(dfSub$avg_log2FC)
    dfOrd = dfSub[order(dfSub$AbsLog, decreasing = T),]
    topMark = head(dfOrd$Genes, topNumb)
    markers = c(markers, topMark)
  }
  unMark = unique(markers)
  return(unMark)
}

topMark = getTopMarkers(df = groupMarkers, topNumb = 5)

RNA.combined.norm$newMonocClust <- factor(RNA.combined.norm$newMonocClust, levels = c("OPC", "Intermideate", "ODC"))

RNA.combined.norm$Group_Cluster = paste(RNA.combined.norm$group, RNA.combined.norm$newMonocClust, sep = "_")

RNA.combined.norm$Group_Cluster <- factor(RNA.combined.norm$Group_Cluster, levels = c("Control_OPC", "Stress_OPC", "Control_Intermideate", "Stress_Intermideate", "Control_ODC", "Stress_ODC"))
  
dotPlot<-DotPlot(object = RNA.combined.norm, features = topMark , scale.max = 100, dot.scale = 16, group.by = 'Group_Cluster')+
  scale_colour_gradient2(low = "blue", mid = "grey", high = "red")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  theme(text = element_text(size = 24))+ # all text elements size
  theme(axis.text = element_text(size = 24))  # axes text size

dotPlot

png(paste0(targDir,'DotPlot_Top5Genes_GroupCluster_', curDate, '.png'), height = 12, width = 24, units = 'in', res = 300)
dotPlot
dev.off()

# AU cell

# load enrichment results
opcOdcCells = rownames(RNA.combined.norm@meta.data)

addEnrich<-function(x){
  load(x)
  AUCmat <- AUCell::getAUC(cells_AUC)
  AucSub = AUCmat[, c(opcOdcCells)]
  RNA.combined.norm[['AUC']] <- CreateAssayObject(data = AucSub)
  DefaultAssay(RNA.combined.norm) <- 'AUC'
  RNA.combined.norm <- ScaleData(RNA.combined.norm, assay = 'AUC', features = rownames(AUCmat))
  return(RNA.combined.norm)
}

# add AUcell data
RNA.combined.norm<-addEnrich(x='./cellsAUC_keggClustProf.RData')
DefaultAssay(RNA.combined.norm)


groupMarkers = read.csv("OPC_ODC/StressVsContr/AUC_keggClustProf_DE/DiffExprMonoc3Clust_AllClust_ContrVsStress_2023-04-10.csv")

topMark = getTopMarkers(df = groupMarkers, topNumb = 5)


groupMarkers = read.csv("OPC_ODC/StressVsContr/RNA_DE/DiffExprMonoc3Clust_AllClust_ContrVsStress_2023-04-10.csv")

RNA.combined.norm$newMonocClust = RNA.combined.norm$MonocClust
RNA.combined.norm$newMonocClust[RNA.combined.norm$newMonocClust == 4] = 3
RNA.combined.norm$newMonocClust[RNA.combined.norm$newMonocClust == 1] = "ODC"
RNA.combined.norm$newMonocClust[RNA.combined.norm$newMonocClust == 2] = "OPC"
RNA.combined.norm$newMonocClust[RNA.combined.norm$newMonocClust == 3] = "Intermideate"
RNA.combined.norm$newMonocClust <- factor(RNA.combined.norm$newMonocClust, levels = c("OPC", "Intermideate", "ODC"))

RNA.combined.norm$Group_Cluster = paste(RNA.combined.norm$group, RNA.combined.norm$newMonocClust, sep = "_")

RNA.combined.norm$Group_Cluster <- factor(RNA.combined.norm$Group_Cluster, levels = c("Control_OPC", "Stress_OPC", "Control_Intermideate", "Stress_Intermideate", "Control_ODC", "Stress_ODC"))

dotPlot<-DotPlot(object = RNA.combined.norm, features = topMark , scale.max = 100, dot.scale = 16, group.by = "Group_Cluster")+
  scale_colour_gradient2(low = "blue", mid = "grey", high = "red")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  theme(text = element_text(size = 16))+ # all text elements size
  theme(axis.text = element_text(size = 16)) + 
  theme(axis.text.x = element_text(size = 14)) +
  scale_x_discrete(limits=rev)

dotPlot

png(paste0(targDir,'DotPlot_Top5Genes_AuCell_KeggClustProf_GroupCluster', curDate, '.png'), height = 12, width = 24, units = 'in', res = 300)
dotPlot
dev.off()

