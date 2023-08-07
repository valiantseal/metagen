library(karyoploteR)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

startR = c(8647705, 8666955)
endR = c(8648095, 8667345)


makePlot = function(curDir, curType, curStart, curEnd, curChr, curAnti, ymax, curCol) {
  regStart = curStart - 5000
  regEnd = curEnd + 5000
  #regStart = format(regStart, big.mark = ",")
  #regEnd = format(regEnd, big.mark = ",")
  curRange = paste0(curChr, ":", regStart, "-", regEnd)
  curRegion <- toGRanges(curRange)
  
  outFile = paste0(curDir, "/", toupper(curAnti), "_", curChr, "-",curStart, "-", curEnd,"_peaks_Edit.png" )
  png(file =  outFile, height = 12, width = 16, units = 'in', res = 300)
  kp <- plotKaryotype(zoom = curRegion, genome = "mm10", cex=2)
  genes.data <- makeGenesDataFromTxDb(TxDb.Mmusculus.UCSC.mm10.knownGene,
                                    karyoplot=kp,
                                    plot.transcripts = T, 
                                    plot.transcripts.structure = T)
  genes.data <- addGeneNames(genes.data)
  genes.data <- mergeTranscripts(genes.data)
  kp <- plotKaryotype(zoom = curRegion, cex=2)
  kpPlotGenes(kp, data=genes.data, r0=0, r1=0.15, gene.name.cex = 2)
  file1 <- paste0(curDir, "/", curType, "1_scaleRc_chr4_edit_sort.bw")
  file2 <- paste0(curDir, "/", curType, "2_scaleRc_chr4_edit_sort.bw")
  
  kp <- kpPlotBigWig(kp, data=file1, r0=0.1, r1=0.5, col = curCol, ymin=0, ymax=ymax)
  computed.ymax <- kp$latest.plot$computed.values$ymax
  kpAxis(kp, ymin=0, ymax=computed.ymax, r0=0.1, r1=0.5, cex = 2)
  kpAddLabels(kp, labels = paste0(toupper(curType),"1"), r0=0.1, r1=0.5, cex=2, label.margin = 0.045)
  
  kp <- kpPlotBigWig(kp, data=file2, r0=0.6, r1=1, col = curCol, ymin=0, ymax=ymax)
  computed.ymax <- kp$latest.plot$computed.values$ymax
  kpAxis(kp, ymin=0, ymax=computed.ymax, r0=0.6, r1=1, cex = 2)
  kpAddLabels(kp, labels = paste0(toupper(curType),"2"), r0=0.6, r1=1, cex=2, label.margin = 0.045)
  
  kpAddBaseNumbers(kp, tick.dist = 2000,
                   add.units = TRUE, cex=2, digits = 2, units = 'b')
  
  title(main = paste0(toupper(curAnti), "_", curChr, "-",curStart, "-", curEnd), cex.main = 2.5)
  dev.off()
  
}

makePlot(curDir = "bowtie2", curType = "c", curStart = startR[1], curEnd = endR[1], curChr = "chr4", curAnti = "Chd7", ymax = 18, curCol = "red")
makePlot(curDir = "bowtie2", curType = "c", curStart = startR[2], curEnd = endR[2], curChr = "chr4", curAnti = "Chd7", ymax = 18, curCol = "red")

makePlot(curDir = "bowtie2", curType = "h", curStart = startR[1], curEnd = endR[1], curChr = "chr4", curAnti = "H3K27AC", ymax = 22, curCol = "blue")
makePlot(curDir = "bowtie2", curType = "h", curStart = startR[2], curEnd = endR[2], curChr = "chr4", curAnti = "H3K27AC", ymax = 22, curCol = "blue")

# find h3k peaks that overlap Chd7

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
  df  = read.table(curFile)
  dfSub = df
  for ( i in 1:nrow(dfSub)) {
    curStart = dfSub$V2[i]
    curEnd = dfSub$V3[i]
    curCheck = checkRange(curStart = curStart, curEnd = curEnd)
    if (any(grepl("TRUE", curCheck))) {
      curSub = dfSub[i,]
      combDat = rbind(combDat, curSub)
    }
    
    
  }
  return(combDat)
}

curDf = compRange("../h3k_chd7_shun.txt")