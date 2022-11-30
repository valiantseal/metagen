library(Seurat)
library(Signac)
library(ggplot2)

setwd("C:/Users/abomb/Projects/AU/wang/SingleCellStress/output")

atacInt<-readRDS('atacIntegrated')

DefaultAssay(atacInt)

#gene.activity<-GeneActivity(atacInt)


atacCells<-rownames(atacInt@meta.data)

source('../../programs/renameClusters.R')

rnaCells<-rownames(RNA.combined.norm@meta.data)

length(atacCells)
length(rnaCells)

length(atacCells[(atacCells%in%rnaCells)])


DefaultAssay(RNA.combined.norm)

RNA.combined.norm[['CellName']]<-colnames(RNA.combined.norm)


#atacAssay<-GetAssay(object = atacInt, assay = 'ATAC')

rnaFilt <- subset(RNA.combined.norm, subset = CellName %in% atacCells )

#rm(RNA.combined.norm, atacInt)
#gc()

#rnaFilt[['ATAC']]<-atacAssay

#DefaultAssay(rnaFilt)<-'ATAC'

#unique(Idents(rnaFilt))

#rm(atacAssay, atacCells, rnaCells)
#gc()

#rnaFilt <- RunUMAP(rnaFilt, reduction = "integrated_lsi", dims = 2:30)

#atacInt <- RunUMAP(atacInt, reduction = "integrated_lsi", dims = 2:30)

rnaFilt <- RunUMAP(rnaFilt, reduction = "lsi", dims = 2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_")


saveRDS(RNA.combined.norm, file = 'integRNADoublFilt')

# annotate atac data based on the seurat guide.
rm(RNA.combined.norm)
gc()

atacInt <- RunUMAP(atacInt, reduction = "integrated_lsi", dims = 2:30)

p1 <- DimPlot(rnaFilt, group.by = "Annotations", label = TRUE) + NoLegend() + ggtitle("RNA")
p2 <- DimPlot(atacInt, group.by = "orig.ident", label = FALSE) + NoLegend() + ggtitle("ATAC")
p1 + p2

atacInt[['RNA']]<-NULL
rnaFilt[['integrated']]<-NULL
gc()

DefaultAssay(atacInt)
DefaultAssay(rnaFilt)

# quantify gene activity
rnaFilt <- FindVariableFeatures(rnaFilt)

gene.activities <- GeneActivity(atacInt, features = VariableFeatures(rnaFilt))

# add gene activities as a new assay
atacInt[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)

# normalize gene activities
DefaultAssay(atacInt) <- "ACTIVITY"
atacInt <- NormalizeData(atacInt)
atacInt <- ScaleData(atacInt, features = rownames(atacInt))


# Identify anchors
transfer.anchors <- FindTransferAnchors(reference = rnaFilt, query = atacInt, features = VariableFeatures(object = rnaFilt),
                                        reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca")

celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = rnaFilt$Annotations,
                                     weight.reduction = atacInt[["integrated_lsi"]], dims = 2:30)

atacInt <- AddMetaData(atacInt, metadata = celltype.predictions)

DimPlot(atacInt, group.by = "predicted.id", label = FALSE)  + ggtitle("ATAC")

atacInt[['seurat_annotations']]<-rnaFilt@meta.data$Annotations

atacInt$annotation_correct <- atacInt$predicted.id == atacInt$seurat_annotations
p1 <- DimPlot(atacInt, group.by = "predicted.id", label = TRUE)  + ggtitle("Predicted annotation")
p2 <- DimPlot(atacInt, group.by = "seurat_annotations", label = TRUE)  + ggtitle("RNA based aanotation")
p1 | p2


saveRDS(atacInt, file = 'atacInt')