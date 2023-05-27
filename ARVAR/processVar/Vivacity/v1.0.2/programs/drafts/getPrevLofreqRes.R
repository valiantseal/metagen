dirList = list.dirs("/home/ubuntu/extraVol/ARVAR/288_Rose/original_files/", recursive = F)

getLofreqRes = function() {
  combDat = data.frame(matrix(nrow = 0, ncol = 0))
  dirList = list.dirs("/home/ubuntu/extraVol/ARVAR/288_Rose/original_files/", recursive = F)
  for (curDir in dirList) {
    fileName = paste0(curDir, "/Lofreq/", basename(curDir), "_lofreq-output.txt")
    try({
      curDf = read.delim(fileName)
      sampleName = gsub("GA-EHC-", "", basename(curDir))
      sampleName = gsub("\\-.*", "", sampleName)
      curDf$Sample = sampleName
      curDf$Samp_Pos_Ref_Alt = paste(curDf$Sample, curDf$POSITION, curDf$REF.NT, curDf$VAR.NT, sep = "__")
      curDf$RawVarFreq = (curDf$FWD.VAR + curDf$REV.VAR) / (curDf$FWD.VAR + curDf$REV.VAR + curDf$FWD.REF + curDf$REV.REF)
      combDat = rbind(combDat, curDf)
      
    })
    
  }
  return(combDat )
}

lofreqRes = getLofreqRes()

q1 = toupper(unique(lofreqRes$Samp_Pos_Ref_Alt))

length(q1[q1%in%allConsSnv])

allConsSnv[!(allConsSnv%in%q1)]
# get the frequencies for the truth set

lofreqPres = lofreqRes[(lofreqRes$Samp_Pos_Ref_Alt%in%allConsSnv),]
lofreqSub = lofreqPres[, c('Samp_Pos_Ref_Alt', "RawVarFreq")]

consSnvDf = data.frame(allConsSnv)
colnames(consSnvDf) = "Samp_Pos_Ref_Alt"

trueSnvs = merge(consSnvDf, lofreqSub, by = "Samp_Pos_Ref_Alt")

length(unique(trueSnvs$Samp_Pos_Ref_Alt))

getMeanRawFreq = function(df) {
  snvs = unique(df$Samp_Pos_Ref_Alt)
  meanDf = data.frame(matrix(nrow = 0, ncol = 0))
  for (snv in snvs) {
    dfSub = df[(df['Samp_Pos_Ref_Alt'] == snv),]
    curMean = data.frame(cbind(snv, mean(dfSub$RawVarFreq)))
    colnames(curMean) = c("Samp_Pos_Ref_Alt", 'RawVarFreq')
    meanDf = rbind(meanDf, curMean)
  }
  return(meanDf)
}

trueSnvsMean = getMeanRawFreq(df = trueSnvs)

trueSnvsMeanFilt = trueSnvsMean[(trueSnvsMean$RawVarFreq >= 0.02) & (trueSnvsMean$RawVarFreq <= 0.98),]

rm(list = ls()[!ls() %in% c("trueSnvsMeanFilt", "allConsSnv", "metaResDf")])
selConsSnv = unique(trueSnvsMeanFilt$Samp_Pos_Ref_Alt)

