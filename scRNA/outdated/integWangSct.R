library(Seurat)
library(SeuratData)
library(patchwork)
library(dplyr)
library(writexl)
library(ggplot2)
library(gridExtra)

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

controlFiltered <- subset(RNA_control, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

RNA_stress[["percent.mt"]] <- PercentageFeatureSet(RNA_stress, pattern = "^mt-")

stressFiltered <- subset(RNA_stress, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#you can use RNA_combined.list as "ifnb.list" in their example
RNA_combined.list = list(control = controlFiltered, stress = stressFiltered)

RNA_combined.list <- lapply(X = RNA_combined.list, FUN = function(x) {
  x <- SCTransform(x)  
  # Don't run ScaleData here
  x <- RunPCA(x, verbose = FALSE)
})

features <- SelectIntegrationFeatures(object.list = RNA_combined.list, nfeatures = 3000)

RNA_combined.list <- PrepSCTIntegration(object.list = RNA_combined.list, anchor.features = features)

RNA.anchors <- FindIntegrationAnchors(object.list = RNA_combined.list, normalization.method = "SCT",
                                      anchor.features = features, reduction = "rpca")
RNA.combined.sct <- IntegrateData(anchorset = RNA.anchors, normalization.method = "SCT")

# remove everything from the memory except of intergrated set
rm(list = ls()[!ls() %in% c("RNA.combined.sct")])

###

#saveRDS(RNA.combined.sct, file = 'intergratedData')
RNA.combined.sct<-readRDS('intergratedData')

RNA.combined.sct <- RunPCA(RNA.combined.sct, verbose = FALSE)
DimPlot(RNA.combined.sct, reduction = "pca")
ElbowPlot(RNA.combined.sct, ndims = 50)

setwd("../output")
ggsave('elbowPlot.jpeg', height = 6, width = 10, units = 'in', dpi = 300)

# clustering
RNA.combined.sct <- FindNeighbors(RNA.combined.sct, dims = 1:40)
RNA.combined.sct <- RunUMAP(RNA.combined.sct, reduction = "pca", dims = 1:40)

RNA.combined.sct <- FindClusters(RNA.combined.sct, resolution = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0))
#DimPlot(RNA.combined.sct, group.by = "integrated_snn_res.0.1", label = T)
resols<-c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)



plotUmpap<-function(x){
  for (i in x){
    current_resol<-paste0("integrated_snn_res.", i)
    Idents(RNA.combined.sct)<-current_resol
    p1 <- DimPlot(RNA.combined.sct, reduction = "umap", group.by = "group", label = F)+ggtitle(current_resol)
    p2 <- DimPlot(RNA.combined.sct, reduction = "umap", label = T, repel=T)+ggtitle(current_resol)
    p3<-grid.arrange(p1,p2, ncol=2, nrow=1)
    ggsave(paste0('umapClustering_',current_resol,".jpeg"), plot=p3, height = 6, width = 10, units = 'in', dpi = 300)
  }
}

#plotUmap(resols)

# expression analysis
DefaultAssay(RNA.combined.sct)
DefaultAssay(RNA.combined.sct) <- "SCT"
DefaultAssay(RNA.combined.sct)
Idents(RNA.combined.sct)<-'integrated_snn_res.0.3'

#q1<-dimnames(RNA.combined.sct@assays$RNA)
#q3<-data.frame(q1[1])


# make lists with gene names
all_select_genes<-c("Rbfox3", "Snap25", "Syt1", "Slc17a6", "Grin2b", "Gad1", "Gad2","Slc1a3", "Gfap", "Aldoc", "Glul", "Rarres2", "Slc6a13", "Plp1", "Vcan", "Hmha1", "Iba-1","Alf-1","CX3CR1", "P2ry12", "Reln", "Ndnf")
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
# merge lists
cell_type_list<-list(all_select_genes,neuron, astrocyte, oligodendrocytes, oligodendrocytes_prec, microglia, HexbC, CA1, CA2, CA3,DG, all_neuron_subtypes)
names_cells<-c('all_select_genes','neuron', 'astrocyte', 'oligodendrocytes', 'oligodendrocytes_prec', 'microglia', 'HexbC', 'CA1', 'CA2', 'CA3','DG', 'all_neuron_subtypes')
names(cell_type_list)<-names_cells

# plot function
exprPlot<-function(x) {
    for (cell_name in x) {
    genes<-cell_type_list[[cell_name]]
    fPlot<-FeaturePlot(RNA.combined.sct, features = genes, min.cutoff = "q9")
    ggsave(paste0('featurePlot_',cell_name,".jpeg"), plot=fPlot, height = 10, width = 16, units = 'in', dpi = 300)
    }
}
# subset
neuron_subtype<-c('CA1', 'CA2', 'CA3','DG')
run_names<-'all_neuron_subtypes'
# run feature plot function

#exprPlot(run_names)


