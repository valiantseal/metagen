library(Seurat)
library(ggplot2)
library(stringr)


setwd("C:/Users/abomb/Projects/AU/wang/SingleCellStress/output")

Markers<-read.csv('../data/Not_found_markers_edit.csv')

#marker1<-strsplit(Markers$X, ",")[[1]]
#marker1_id<-Markers[1,1]
#marker1_comb<-c(marker1_id, marker1)
#marker1Filtr<-gsub(" ", "", marker1_comb)

editMarkers<-character()

for (i in 1:nrow(Markers)){
  marker1<-strsplit(Markers$X, ",")[[i]]
  marker1_id<-Markers[i,1]
  marker1_comb<-c(marker1_id, marker1)
  marker1Filtr<-gsub(" ", "", marker1_comb)
  editMarkers<-c(editMarkers, marker1Filtr)
  editMarkUn<-unique(editMarkers)
}



RNA.combined.norm<-readRDS('intergratedRNA')


DefaultAssay(RNA.combined.norm) <- "RNA"

DefaultAssay(RNA.combined.norm) <- "RNA"
DefaultAssay(RNA.combined.norm)
Idents(RNA.combined.norm)<-'integrated_snn_res.0.3'

dimNames<-RNA.combined.norm@assays$RNA@data@Dimnames
all_genes<-dimNames[[1]]

markersPresent<-editMarkUn[(editMarkUn%in%all_genes)]



FeaturePlot(RNA.combined.norm, features = markersPresent, min.cutoff = "q9")
ggsave("Arh_Cx_NewFoundMarkers.jpeg", height = 10, width = 16, units = 'in', dpi = 300)

markersTitle<-str_to_title(tolower(editMarkUn))
markersTitlePresnet<-markersTitle[(markersTitle%in%all_genes)]

markersTitlePresnet1<-markersTitlePresnet[1:8]

markersTitlePresnet2<-markersTitlePresnet[8:15]

markersTitlePresnet3<-markersTitlePresnet[16:23]

markersTitlePresnet4<-markersTitlePresnet[24:32]

FeaturePlot(RNA.combined.norm, features = markersTitlePresnet4, min.cutoff = "q9")
ggsave("TitlePresnet4_NewFoundMarkers.jpeg", height = 10, width = 16, units = 'in', dpi = 300)


# more markers
group1<-c('Trps1', 'Calcrl', 'Dcn', 'Meis2', 'Nectin3', 'St8sia6', 'Fn1', 'Col6a1')
group2<-str_to_title(tolower(c('CSPG4','RNF5','EOMES','SOX2','SOX1','PAX3','PAX6','OTX2','CNTNAP1','ASCL1','SMARCA4','MSI1','MSI2', 'SOX9'))) # 'NES' expression is 0 for all cells
group2a<-group2[1:6]
group2b<-group2[7:14]

FeaturePlot(RNA.combined.norm, features = group2, min.cutoff = "q9")
ggsave("Additional_markers_group2.jpeg", height = 10, width = 16, units = 'in', dpi = 300)