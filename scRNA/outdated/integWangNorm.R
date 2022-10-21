library(Seurat)
library(SeuratData)
library(patchwork)
library(dplyr)
library(writexl)
library(ggplot2)
library(gridExtra)
library(clustree)
library(openxlsx)

#please set the working directory after cleaning the memory
setwd("C:/Users/abomb/Projects/AU/wang/SingleCellStress/data")
#setwd('C:/Users/shyan/Box/Wang_and_Jiao_sequencing_projects/Project1_Wang_single_cell_sequencing/Andrei/')

#please set the working directory after cleaning the memory
setwd("C:/Users/abomb/Projects/AU/wang/SingleCellStress/data")

#here the 10X function is their updated function to read the data, it can read RNA-Seq and ATAC-Seq data together - you need to provide the path to "filtered_feature_bc_matrix" folder in the output folder of Cellranger)
control = Read10X(data.dir = "../data/control")

stress = Read10X(data.dir = "../data/stress")

#create SeuratObject for RNA-Seq data
RNA_control = CreateSeuratObject(counts = control$`Gene Expression`)
RNA_control$group = "Control" #create a group label to identify the two groups

RNA_stress = CreateSeuratObject(counts = stress$`Gene Expression`)
RNA_stress$group = "Stress" #create a group label to identify the two groups

# quality control 
RNA_control[["percent.mt"]] <- PercentageFeatureSet(RNA_control, pattern = "^mt-")

median(RNA_control$nFeature_RNA)
sd(RNA_control$nFeature_RNA)

VlnPlot(RNA_control, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

controlFiltered <- subset(RNA_control, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

RNA_stress[["percent.mt"]] <- PercentageFeatureSet(RNA_stress, pattern = "^mt-")

median(RNA_stress$nFeature_RNA)

VlnPlot(RNA_stress, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

stressFiltered <- subset(RNA_stress, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#you can use RNA_combined.list as "ifnb.list" in their example
RNA_combined.list = list(control = controlFiltered, stress = stressFiltered)




# normalize and identify variable features for each dataset independently
RNA_combined.list <- lapply(X = RNA_combined.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = RNA_combined.list)
immune.anchors <- FindIntegrationAnchors(object.list = RNA_combined.list, anchor.features = features)
RNA.combined.norm <- IntegrateData(anchorset = immune.anchors)
rm(list = ls()[!ls() %in% c("RNA.combined.norm")])

DefaultAssay(RNA.combined.norm)
RNA.combined.norm  <- ScaleData(RNA.combined.norm , verbose = FALSE)
RNA.combined.norm <- RunPCA(RNA.combined.norm, verbose = FALSE)

#saveRDS(RNA.combined.norm, file = 'intergratedRNA')
#RNA.combined.norm<-readRDS('intergratedRNA')

setwd("../output")

# clustering
DimPlot(RNA.combined.norm, reduction = "pca")
ElbowPlot(RNA.combined.norm, ndims = 50)

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

minPCA<-fgetNPCs(obj=RNA.combined.norm, MIN_CUM_SD=90)

#RNA.combined.norm <- JackStraw(RNA.combined.norm, num.replicate = 100, dims = 50)
#RNA.combined.norm <- ScoreJackStraw(RNA.combined.norm, dims = 1:50)
#JackStrawPlot(RNA.combined.norm, dims = 1:50)

RNA.combined.norm <- FindNeighbors(RNA.combined.norm, dims = 1:40)
RNA.combined.norm <- RunUMAP(RNA.combined.norm, reduction = "pca", dims = 1:40)

RNA.combined.norm <- FindClusters(RNA.combined.norm, resolution = c(0, 0.05, 0.08,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0))
#DimPlot(RNA.combined.norm, group.by = "integrated_snn_res.0.1", label = T)
resols<-c(0,0.05,0.08,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)

# terst clusters
clustree <- FindClusters(RNA.combined.norm , resolution = resols)
p <- clustree(clustree)
ggsave("Clustree.jpeg", width = 8, height = 10)

plotUmap<-function(x){
  for (i in x){
    current_resol<-paste0("integrated_snn_res.", i)
    Idents(RNA.combined.norm)<-current_resol
    p1 <- DimPlot(RNA.combined.norm, reduction = "umap", group.by = "group", label = F)+ggtitle(current_resol)
    p2 <- DimPlot(RNA.combined.norm, reduction = "umap", label = T, repel=T)+ggtitle(current_resol)
    p3<-grid.arrange(p1,p2, ncol=2, nrow=1)
    ggsave(paste0('umapClustering_',current_resol,".jpeg"), plot=p3, height = 6, width = 10, units = 'in', dpi = 300)
  }
}

plotUmap(resols)

# differential expression
DefaultAssay(RNA.combined.norm) <- "RNA"
DefaultAssay(RNA.combined.norm)
Idents(RNA.combined.norm)<-'integrated_snn_res.0.3'

# check for normalization
q1<-RNA.combined.norm@assays$RNA@counts
q2<-RNA.combined.norm@assays$RNA@data
identical(q1,q2)

RNA.combined.norm <- ScaleData(RNA.combined.norm, verbose = FALSE)
# save with all plots and modifications
saveRDS(RNA.combined.norm, file = 'intergratedRNA')
RNA.combined.norm<-readRDS('intergratedRNA')
###

# make lists with gene names
all_select_genes<-c("Rbfox3", "Snap25", "Syt1", "Slc17a6", "Grin2b", "Gad1", "Gad2","Slc1a3", "Gfap", "Aldoc", 
                    "Glul", "Rarres2", "Slc6a13", "Plp1", "Vcan", "Hmha1", "Iba-1","Alf-1","CX3CR1", "P2ry12", "Reln", "Ndnf")
neuron<-c("Rbfox3", "Snap25", "Syt1", "Slc17a6", "Grin2b", "Gad1", "Gad2")
astrocyte<-c("Slc1a3", "Gfap", "Aldoc", "Glul", "Rarres2", "Slc6a13")
oligodendrocytes<-c("Plp1")
oligodendrocytes_prec<-c("Vcan")
microglia<-c("Hmha1", "Iba-1","Alf-1","CX3CR1", "P2ry12")
HexbC<- c("Reln", "Ndnf")
CA1<-c("Mpped1", "Satb2")
CA2<-c("Map3k15", "Gram4")
CA3<-c("Cdh24", "Npy2r")
DG<-c("Prox1", "Glis3")
all_neuron_subtypes<-c("Reln", "Ndnf", "Mpped1", "Satb2", "Map3k15", "Gram4","Cdh24", "Npy2r",  "Prox1", "Glis3")
additional_markers<-c("Cx3cr1", "Ephb1", "Vwc2I", "Csf2rb2", "Fibcd1", "Wfs1","DCN",  "Slc4a4", "Plk5" , "Cntn6", "Amigo2")
addMarkSet2<-c("Kcnh5", "Bok", "CSmd3","PCDH9", "Akt2", "Acsl3", "Trf", "Ptgds", "Hexb", "Kcnc1", "Npy1r", "Calb2", "Calb1")
addMarkSet3<-c("Olig2", "Gramd4", "Vwc2l")
# merge lists
cell_type_list<-list(neuron, astrocyte, oligodendrocytes, oligodendrocytes_prec, microglia, HexbC, CA1, CA2, CA3,DG, all_neuron_subtypes, 
                     additional_markers, addMarkSet2, addMarkSet3)
names_cells<-c('neuron', 'astrocyte', 'oligodendrocytes', 'oligodendrocytes_prec', 'microglia', 'HexbC', 'CA1', 'CA2', 'CA3','DG', 'all_neuron_subtypes', 
               'additional_markers', 'addMarkSet2', 'addMarkSet3')
names(cell_type_list)<-names_cells

# plot function
exprPlot<-function(x) {
  for (cell_name in x) {
    genes<-cell_type_list[[cell_name]]
    fPlot<-FeaturePlot(RNA.combined.norm, features = genes, min.cutoff = "q9")
    ggsave(paste0('featurePlot_',cell_name,".jpeg"), plot=fPlot, height = 10, width = 16, units = 'in', dpi = 300)
  }
}

exprPlot(names_cells)

FeaturePlot(RNA.combined.norm, features = "lyd", min.cutoff = "q9")
ggsave('./output/featurePlot_addMarkSet3.jpeg', height = 10, width = 16, units = 'in', dpi = 300)


# cluster analysis
metadata<-data.frame(RNA.combined.norm@meta.data)

cellsPerClust<-data.frame(table(metadata$integrated_snn_res.0.3, metadata$group))
colnames(cellsPerClust)<-c('Cluster', 'Group', 'Cells_number')
write.xlsx(cellsPerClust, './output/Cells_per_cluster_group_res0.3.xlsx')


# more markers
pyr_markers<-read.csv('../data/pyr_markers.csv')

pyr_set1<-pyr_markers[1:10,]
pyr_set2<-pyr_markers[11:20,]
pyr_set3<-pyr_markers[21:30,]


pyr_list<-list(pyr_set1, pyr_set2, pyr_set3)
names(pyr_list)<-c('pyr_set1', 'pyr_set2', 'pyr_set3')
pyr_names<-c('pyr_set1', 'pyr_set2', 'pyr_set3')

# plot function
exprPlot<-function(x, y) {
  for (cell_name in x) {
    genes<-y[[cell_name]]
    fPlot<-FeaturePlot(RNA.combined.norm, features = genes, min.cutoff = "q9")
    ggsave(paste0('featurePlot_',cell_name,".jpeg"), plot=fPlot, height = 10, width = 16, units = 'in', dpi = 300)
  }
}

exprPlot(pyr_names, pyr_list)
