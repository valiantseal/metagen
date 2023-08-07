library(ggplot2)

library(karyoploteR)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)


curRegion <- toGRanges("chr4:8848930-8854067")

kp <- plotKaryotype(zoom = curRegion, genome = "mm10", cex=2)

genes.data <- makeGenesDataFromTxDb(TxDb.Mmusculus.UCSC.mm10.knownGene,
                                    karyoplot=kp,
                                    plot.transcripts = TRUE, 
                                    plot.transcripts.structure = TRUE)




kp <- plotKaryotype(zoom = curRegion, cex=2, genome = "mm10")
kpAddBaseNumbers(kp, tick.dist = 1000, minor.tick.dist = 100,
                 add.units = TRUE, cex=1.3, digits = 6)
kpPlotGenes(kp, data=genes.data, r0=0, r1=0.15, gene.name.cex = 2)

c1 <- "c1_figs/c1_rpkmNorm_bin1.bw"
kp <- kpPlotBigWig(kp, data=c1, r0=0.35, r1=0.65, col = "blue", ymin=0, ymax=1000)
computed.ymax <- kp$latest.plot$computed.values$ymax
kpAxis(kp, ymin=0, ymax=computed.ymax, r0=0.35, r1=0.65)
kpAddLabels(kp, labels = "C1", r0=0.22, r1=0.3, cex=2)


i1 <- "c1_figs/i1_rpkmNorm_bin1.bw"
kp <- kpPlotBigWig(kp, data=i1, ymax=700, r0=0.7, r1=1, ymin = 0, col = "black")
computed.ymax <- kp$latest.plot$computed.values$ymax
kpAxis(kp, ymin=0, ymax=computed.ymax, r0=0.7, r1=1)
kpAddLabels(kp, labels = "I1", r0=1.3, r1=1, cex=1.6, label.margin = 0.035)

print(kp)

##
png("c1_i1_karyopTr.png", height = 16, width = 24, units = "in", res = 300)
curRegion <- toGRanges("chr4:8848930-8854067")

kp <- plotKaryotype(zoom = curRegion, genome = "mm10", cex=2)

genes.data <- makeGenesDataFromTxDb(TxDb.Mmusculus.UCSC.mm10.knownGene,
                                    karyoplot=kp,
                                    plot.transcripts = TRUE, 
                                    plot.transcripts.structure = TRUE)

kp <- plotKaryotype(zoom = curRegion, cex=2, genome = "mm10")
kpAddBaseNumbers(kp, tick.dist = 1000, minor.tick.dist = 100,
                 add.units = TRUE, cex=1.3, digits = 6)

c1 <- "c1_figs/c1_rpkmNorm_bin1.bw"
kp <- kpPlotBigWig(kp, data=c1, ymax="visible.region", r0=0.35, r1=0.65, col = "blue")
computed.ymax <- kp$latest.plot$computed.values$ymax
kpAxis(kp, ymin=0, ymax=computed.ymax, r0=0.35, r1=0.65)
kpAddLabels(kp, labels = "C1", r0=0.22, r1=0.3, cex=2)


i1 <- "c1_figs/i1_rpkmNorm_bin1.bw"
kp <- kpPlotBigWig(kp, data=i1, ymax="visible.region", r0=0.7, r1=1)
computed.ymax <- kp$latest.plot$computed.values$ymax
kpAxis(kp, ymin=0, ymax=computed.ymax, r0=0.7, r1=1)
kpAddLabels(kp, labels = "I1", r0=1.3, r1=1, cex=1.6, label.margin = 0.035)

print(kp)

dev.off()