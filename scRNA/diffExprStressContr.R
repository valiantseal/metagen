# should follow after renameClust.R

library(Seurat)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(clustree)

curDate<-Sys.Date()

DefaultAssay(RNA.combined.norm) <- "RNA"

DefaultAssay(RNA.combined.norm)
levels(RNA.combined.norm)

clusters<-c('CA1', 'CA2', 'CA3', 'DG', 'GABA', 'ODC', 'OPC', 'MG','CR')

RNA.combined.norm[["Annotations"]] <- Idents(object = RNA.combined.norm)

# all markers 


cellType<-subset(x = RNA.combined.norm, subset = Annotations == "CA1")

Idents(cellType)<-cellType$group

allMarkers <- FindMarkers(cellType, ident.1 = "Control", ident.2 = "Stress", only.pos = F, min.pct = 0.1, logfc.threshold = 0, test.use = "MAST")

allMarkers$Genes<-rownames(allMarkers)

allMarkers<-tibble::add_column(allMarkers, Cell_Type='CA1', .before = 1)



findMarkersGr<-function(dat, clust, pos){
  combMarkers<-data.frame(matrix(ncol=0, nrow=1))
  for ( i in clust){
    cellType<-subset(x = dat, subset = Annotations == i)
    Idents(cellType)<-cellType$group
    allMarkers <- FindMarkers(cellType , only.pos = pos, min.pct = 0.1, logfc.threshold = 0, test.use = "MAST", ident.1 = "Control", ident.2 = "Stress")
    allMarkers$Genes<-rownames(allMarkers)
    write.csv(allMarkers, paste0('DiffExpr_', i, '_ContrVsStress_', curDate, '.csv'), row.names = F)
    allMarkers<-tibble::add_column(allMarkers, Cell_Type=i, .before = 1)
    combMarkers<-rbind(combMarkers, allMarkers)
  }
  return(combMarkers)
}

groupMarkers<-findMarkersGr(dat=RNA.combined.norm, clust=clusters, pos=F)
write.csv(groupMarkers, paste0('allMarkRenClust0.3ResGroups_',curDate, '.csv'), row.names = F)

#posGroupMarkers<-findMarkersGr(dat=RNA.combined.norm, clust=clusters, pos=T)
#write.csv(posGroupMarkers, paste0('positiveMarkRenClust0.3ResGroups_',curDate, '.csv'), row.names = F)