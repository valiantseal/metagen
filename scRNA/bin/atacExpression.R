library(Seurat)
library(Signac)
library(ggplot2)

setwd("C:/Users/abomb/Projects/AU/wang/SingleCellStress/output")

atacInt<-readRDS('atacIntegrated')

DefaultAssay(atacInt)

gene.activity<-GeneActivity(atacInt)


atacCells<-rownames(atacInt@meta.data)

source('../../programs/renameClusters.R')

rnaCells<-rownames(RNA.combined.norm@meta.data)

length(atacCells)
length(rnaCells)

length(atacCells[(atacCells%in%rnaCells),])