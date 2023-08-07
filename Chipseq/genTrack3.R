library(karyoploteR)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)


checkRange = function(curStart, curEnd) {
  startR = c(8647705, 8666955)
  endR = c(8648095, 8667345)
  combDat = character()
  
  for ( i in 1:length(startR) ) {
    x1 = startR[i]
    x2 = endR[i]
    
    checkVar = ((x1 <= curEnd) && (curStart <= x2))
    checkVar = paste0(i, "_", checkVar)
    combDat = c(combDat, checkVar)
  }
  return(combDat)
}


compRange = function(curFile) {
  combCheck = character()
  combDat = data.frame(matrix(nrow = 0, ncol = 0))
  df  = read.delim(curFile)
  dfSub = df[(df$Gene.Name == "Chd7"),]
  for ( i in 1:nrow(dfSub)) {
    curStart = dfSub$Star[i]
    curEnd = dfSub$End[i]
    curCheck = checkRange(curStart = curStart, curEnd = curEnd)
    if (any(grepl("TRUE", curCheck))) {
      curSub = dfSub[i,]
      combDat = rbind(combDat, curSub)
    }
    
    
  }
  return(combDat)
}

macs2_poolAll = compRange(curFile = "/home/flyhunter/Kai/Chipseq/pnas/process_custom/macs2_poolAll/KJ-C.annotpeaks")

# make plots
targDir = "bamComp/"
dir.create(targDir)




# for Regions overlapping with David's peaks
for ( i in 1:nrow(macs2_poolAll)) {
  regStart = macs2_poolAll$Start[i] - 10000
  regStart = format(regStart, big.mark = ",")
  regEnd = macs2_poolAll$End[i] + 10000
  regEnd = format(regEnd, big.mark = ",")
  chr = macs2_poolAll$Chr[i]
  curRange = paste0(chr, ":", regStart, "-", regEnd)
  curRegion <- toGRanges(curRange)
  
  png(file = paste0("bamComp/C_merge_adjInp_log2_Chd7_peak", i, ".png"), height = 12, width = 16, units = 'in', res = 300)
  
  kp <- plotKaryotype(zoom = curRegion, genome = "mm10", cex=2)
  genes.data <- makeGenesDataFromTxDb(TxDb.Mmusculus.UCSC.mm10.knownGene,
                                      karyoplot=kp,
                                      plot.transcripts = T, 
                                      plot.transcripts.structure = T)
  genes.data <- addGeneNames(genes.data)
  genes.data <- mergeTranscripts(genes.data)
  
  kp <- plotKaryotype(zoom = curRegion, cex=2, genome = "mm10")
  kpPlotGenes(kp, data=genes.data, r0=0, r1=0.15, gene.name.cex = 2)
  
  c1 <- "scale_input/C_merged_scaleRc_log2.bw"
  kp <- kpPlotBigWig(kp, data=c1, r0=0.3, r1=1, col = "black", ymin=0, ymax="visible.region")
  computed.ymax <- kp$latest.plot$computed.values$ymax
  kpAxis(kp, ymin=0, ymax=computed.ymax, r0=0.3, r1=1, cex = 1.5)
  
  # kpAddBaseNumbers(kp, tick.dist = 10000, minor.tick.dist = 2000,
  #                  add.units = TRUE, cex=1.3, digits = 6)
  
  kpAddBaseNumbers(kp, tick.dist = 5000,
                   add.units = TRUE, cex=1.3, digits = 2, units = 'b')
  
  
  title(main = paste0("Chd7 ", curRange), cex.main = 2.5)
  dev.off()
  
}

# Whole Chd7 Regions
regStart = 8691365 - 10000
regStart = format(regStart, big.mark = ",")
regEnd = 8867659 + 10000
regEnd = format(regEnd, big.mark = ",")
chr = "chr4"
curRange = paste0(chr, ":", regStart, "-", regEnd)
curRegion <- toGRanges(curRange)
#curRegion <- toGRanges("chr4:8,691,365-8,867,659")


png(file = paste0("bamComp/C_Merged_adjInp_log2_Whole_Chd7.png"), height = 12, width = 16, units = 'in', res = 300)

kp <- plotKaryotype(zoom = curRegion, genome = "mm10", cex=2)
genes.data <- makeGenesDataFromTxDb(TxDb.Mmusculus.UCSC.mm10.knownGene,
                                    karyoplot=kp,
                                    plot.transcripts = T, 
                                    plot.transcripts.structure = T)
genes.data <- addGeneNames(genes.data)
genes.data <- mergeTranscripts(genes.data)

kp <- plotKaryotype(zoom = curRegion, cex=2)
kpPlotGenes(kp, data=genes.data, r0=0, r1=0.15, gene.name.cex = 2)

c1 <- "scale_input/C_merged_scaleRc_log2.bw"
kp <- kpPlotBigWig(kp, data=c1, r0=0.3, r1=1, col = "black", ymin=0, ymax="visible.region")
computed.ymax <- kp$latest.plot$computed.values$ymax
kpAxis(kp, ymin=0, ymax=computed.ymax, r0=0.3, r1=1, cex = 1.5)

# kpAddBaseNumbers(kp, tick.dist = 10000, minor.tick.dist = 2000,
#                  add.units = TRUE, cex=1.3, digits = 6)

kpAddBaseNumbers(kp, tick.dist = 20000,
                 add.units = TRUE, cex=1.3, digits = 2, units = 'Kb')


title(main = paste0("Chd7"), cex.main = 2.5)
dev.off()
