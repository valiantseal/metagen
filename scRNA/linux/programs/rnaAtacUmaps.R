library(ArchR)
library(stringr)

targDir = './RnaAtacPlots/'

dir.create(targDir)

# get RNA seq UMAP
source('../programs/renameClusters.R')

curDate = Sys.Date()

clusters<-c('SUB', 'CA1', 'CA2', 'CA3', 'DG', 'GABA', 'C-R', 'OPC', 'ODC', 'MG')

df<-data.frame(RNA.combined.norm@reductions$umap@cell.embeddings)
df$Cell_id<-row.names(df)

metadata<-data.frame(RNA.combined.norm@meta.data)
metadata$Cell_id<-row.names(metadata)

identical(metadata$Cell_id, df$Cell_id)

# atac

projFilt <- readRDS('atacArchRnaFilt')
archUmap <- getEmbedding(ArchRProj = projFilt, embedding = "UMAP", returnDF = TRUE)
archUmap$Cell_id = rownames(archUmap)
identical(archUmap$Cell_id, projFilt$cellNames)
archUmap$Annotations = projFilt$Annotations
archUmap$Group = str_to_title(projFilt$Group)
colnames(archUmap)[1:2] = c('UMAP_1', 'UMAP_2')

levels(x = archUmap) <- c('SUB', 'CA1', 'CA2', 'CA3', 'DG', 'GABA', 'C-R', 'OPC', 'ODC', 'MG')
#archUmap$Annotations <- factor(archUmap$Annotations, levels=c('SUB', 'CA1', 'CA2', 'CA3', 'DG', 'GABA', 'C-R', 'OPC', 'ODC', 'MG'))
# plots
plot1 = ggplot(aes(x = UMAP_1, y = UMAP_2, color = Annotations), data = archUmap )+
  geom_point() +
  theme_classic() +
  theme(text = element_text(size = 20)) +
  guides(color = guide_legend(override.aes = list(size = 9))) +
  guides(shape = guide_legend(override.aes = list(size = 7)))
  

plot1

ggsave(plot = plot1, file = paste0(targDir, 'AtacUmapTileMat_allClusersAnnotations_', curDate, '.png'), width = 12, height = 12, dpi = 300, units = 'in')
