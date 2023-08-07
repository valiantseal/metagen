

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
  df  = read.delim(curFile)
  dfSub = df[(df$Gene.Name == "Chd7"),]
  for ( i in 1:nrow(dfSub)) {
    curStart = dfSub$Star[i]
    curEnd = dfSub$End[i]
    curCheck = checkRange(curStart = curStart, curEnd = curEnd)
    combCheck = c(combCheck, curCheck)
    
  }
  print(combCheck)
}

macs2_poolAll = compRange(curFile = "/home/flyhunter/Kai/Chipseq/pnas/process_custom/macs2_poolAll/KJ-C.annotpeaks")
unique(macs2_poolAll)
grep("TRUE", macs2_poolAll)

macs2_poolContr = compRange(curFile = "/home/flyhunter/Kai/Chipseq/pnas/process_custom/macs2_poolContr/KJ-C_merge.annotpeaks")
unique(macs2_poolContr)
grep("TRUE", macs2_poolContr)



chip_mergeContr = compRange(curFile = "/home/flyhunter/Kai/Chipseq/pnas/output_merge_control/bwa/mergedLibrary/macs2/narrowPeak/consensus/CHD7/CHD7.consensus_peaks.annotatePeaks.txt")
unique(chip_mergeContr)
grep("TRUE", chip_mergeContr)

##
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
macs2_poolAll$Start
macs2_poolAll$End


macs2_poolContr = compRange(curFile = "/home/flyhunter/Kai/Chipseq/pnas/process_custom/macs2_poolContr/KJ-C_merge.annotpeaks")
macs2_poolContr$Start
macs2_poolContr$End

macs2_fdr = compRange(curFile = "/home/flyhunter/Kai/Chipseq/pnas/process_custom/macs2_poolAll_fdr0.01/KJ-C.annotpeaks")
macs2_fdr$Start
macs2_fdr$End

macs2_fdr_broad = compRange(curFile = "/home/flyhunter/Kai/Chipseq/pnas/process_custom/macs2_poolAll_fdr0.01_broad/KJ-C.annotpeaks")

macs2_broad = compRange(curFile = "/home/flyhunter/Kai/Chipseq/pnas/process_custom/macs2_poolAll_broad/KJ-C.annotpeaks")

macs2_h = compRange(curFile = "/home/flyhunter/Kai/Chipseq/pnas/process_custom/macs2_poolAll_fdr0.01/KJ-H.annotpeaks")
