library(ggplot2)

library(karyoploteR)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)


curRegion <- toGRanges("chr4:8848930-8854067")

kp <- plotKaryotype(zoom = curRegion, genome = "mm10", cex=2)


genes.data <- makeGenesDataFromTxDb(TxDb.Mmusculus.UCSC.mm10.knownGene,
                                    karyoplot=kp,
                                    plot.transcripts = T, 
                                    plot.transcripts.structure = T)

genes.data <- addGeneNames(genes.data)
genes.data <- mergeTranscripts(genes.data)

kp <- plotKaryotype(zoom = curRegion, cex=2)
kpPlotGenes(kp, data=genes.data, r0=0, r1=0.15, gene.name.cex = 2)

c1 <- "c1_figs/c1_rpkmNorm_bin1.bw"
kp <- kpPlotBigWig(kp, data=c1, r0=0.3, r1=1, col = "blue", ymin=0, ymax=700)

kpAxis(kp, ymin=0, ymax=computed.ymax, r0=0.3, r1=1, cex = 1.5)

# kpAddBaseNumbers(kp, tick.dist = 10000, minor.tick.dist = 2000,
#                  add.units = TRUE, cex=1.3, digits = 6)

kpAddBaseNumbers(kp, tick.dist = 1000, minor.tick.dist = 100,
                 add.units = TRUE, cex=1.3, digits = 6)



 title(main = "Scatter Plot Example")