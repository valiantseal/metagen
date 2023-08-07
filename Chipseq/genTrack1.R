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

# c1 <- "c1_figs/c1_rpkmNorm_bin1.bw"
# kp <- kpPlotBigWig(kp, data=c1, r0=0.1, r1=0.5, col = "blue", ymin=0, ymax=700)
# computed.ymax <- kp$latest.plot$computed.values$ymax
# kpAxis(kp, ymin=0, ymax=computed.ymax, r0=0.1, r1=0.5, cex = 1.5)
# kpAddLabels(kp, labels = "C1", r0=0.1, r1=0.5, cex=2, label.margin = 0.05)
# 
# 
# i1 <- "c1_figs/i1_rpkmNorm_bin1.bw"
# kp <- kpPlotBigWig(kp, data=i1, ymax=700, r0=0.6, r1=1, ymin = 0, col = "black")
# computed.ymax <- kp$latest.plot$computed.values$ymax
# kpAxis(kp, ymin=0, ymax=computed.ymax, r0=0.6, r1=1, cex = 1.5)
# kpAddLabels(kp, labels = "I1", r0=0.6, r1=1, cex=2, label.margin = 0.05)

# try to do peaks merged



c1 <- "c1_figs/c1_rpkmNorm_bin1.bw"
kp <- kpPlotBigWig(kp, data=c1, r0=0.1, r1=1, col = "blue", ymin=0, ymax=700)
computed.ymax <- kp$latest.plot$computed.values$ymax
kpAxis(kp, ymin=0, ymax=computed.ymax, r0=0.1, r1=1, cex = 1.5)
#kpAddLabels(kp, labels = "C1", r0=0.1, r1=1, cex=2, label.margin = 0.05)

i1 <- "c1_figs/i1_rpkmNorm_bin1.bw"
kp <- kpPlotBigWig(kp, data=i1, ymax=700, r0=0.1, r1=1, ymin = 0, col = "black")
computed.ymax <- kp$latest.plot$computed.values$ymax


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
kp <- kpPlotBigWig(kp, data=c1, r0=0.1, r1=1, col = "blue", ymin=0, ymax=700)
computed.ymax <- kp$latest.plot$computed.values$ymax
kpAxis(kp, ymin=0, ymax=computed.ymax, r0=0.1, r1=1, cex = 1.5)
#kpAddLabels(kp, labels = "C1", r0=0.1, r1=1, cex=2, label.margin = 0.05)

i1 <- "c1_figs/i1_rpkmNorm_bin1.bw"
kp <- kpPlotBigWig(kp, data=i1, ymax=700, r0=0.1, r1=1, ymin = 0, col = "black")
dev.off()