library(Seurat)
library(Signac)
library(ggplot2)
library(gridExtra)
library(EnsDb.Mmusculus.v79)

curDate<-Sys.Date()

setwd("/home/flyhunter/Wang/output")

atacInt<-readRDS('atacIntegrated_macs2')

DefaultAssay(atacInt)



clusters<-c('SUB', 'CA1', 'CA2', 'CA3', 'DG', 'GABA', 'C-R', 'OPC', 'ODC', 'MG')

DimPlot(atacInt)

#
# annotation
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)

seqlevelsStyle(annotations) <- 'UCSC'

Annotation(atacInt) <- annotations

# all markers
allMarkers <- FindAllMarkers(atacInt , only.pos = T, min.pct = 0.05, test.use = "LR", 
                          latent.vars = 'atac_peak_region_fragments')
closGenes<-ClosestFeature(atacInt, regions = rownames(allMarkers))

allMarkEdit<-cbind(allMarkers, closGenes)


dir.create('atacIntDE')

targDir<-'atacIntDE/allMarkers/'

dir.create(targDir)

write.csv(allMarkEdit, file = paste0(targDir, 'clustSpecificMarkers.csv'), row.names = F)

for ( i in clusters ){
  dfSel<-allMarkEdit[(allMarkEdit$cluster==i),]
  write.csv(dfSel, file = paste0(targDir, i, '_clustSpecificMarkers.csv'), row.names = F)
}

# function to claculate differential gene expression
findMarkersGr<-function(dat, clust, pos, latVar){
  combMarkers<-data.frame(matrix(ncol=0, nrow=0))
  for ( i in clust){
    # subset cell type from all data
    cellType<-subset(x = dat, subset = Annotations == i)
    # change identity from cell type to group
    Idents(cellType)<-cellType$dataset
    # calculate gene expression
    allMarkers <- FindMarkers(cellType , only.pos = pos, min.pct = 0.05, test.use = "LR", 
                              latent.vars = latVar, 
                              ident.1 = "Stress", ident.2 = "Control")
    closGenes<-ClosestFeature(atacInt, regions = rownames(allMarkers))
    
    allMarkEdit<-cbind(allMarkers, closGenes)
    allMarkers<-tibble::add_column(allMarkers, Cell_Type=i, .before = 1)
    write.csv(allMarkEdit, paste0('./atacIntDE/atacInt_DE_', i, '_ContrVsStress_nolatvar_', curDate, '.csv'), row.names = F)
    combMarkers<-rbind(combMarkers, allMarkEdit)
  }
  return(combMarkers)
}

groupMarkers<-findMarkersGr(dat=atacInt, clust=clusters, pos=F, latVar = NULL) # atac_peak_region_fragments

write.csv(groupMarkers, paste0('./atacIntDE/allDe_atacInt_StrVsContr_nolatvar_',curDate, '.csv'), row.names = F)

# save annotated 
saveRDS(atacInt, file = 'atacIntegrated_macs2')

CoveragePlot(atacInt, region = allMarkEdit$query_region[2],
             extend.upstream = 10000,
             extend.downstream = 5000,
             group.by = "dataset")
