# get vivacity data
getResBind = function(path, freqCol, minFreq, maxFreq) {
  combTab = data.frame(matrix(nrow = 0 , ncol = 0))
  filesList = list.files(path,)
  for ( i in filesList) {
    try({
      inFile = paste0(path, i, "/filtered.csv")
      df = read.csv(inFile)
      df$RawVarFreq = (df$FWD.VAR + df$REV.VAR) / (df$FWD.VAR + df$REV.VAR + df$FWD.REF + df$REV.REF)
      dfFilt = df[(df[freqCol] >= minFreq) & (df[freqCol] <= maxFreq),]
      dfFilt$ExactSamp = i
      sampleList = strsplit(i, '-')
      sampleName =  sampleList[[1]][1:3]
      sampleName = paste(sampleName, collapse = "-")
      dfFilt$Sample = sampleName
      dfFilt$Samp_Pos_Ref_Alt = paste(dfFilt$Sample, dfFilt$POSITION, dfFilt$REF.NT, dfFilt$VAR.NT, sep = "__")
      combTab = rbind(combTab, dfFilt)
    })

  }
  return(combTab)
}


# calculate strand bias ratio for variant 
addVarSB = function(metaResDf) {
  metaResDf$Var_SB = metaResDf$FWD.VAR / metaResDf$REV.VAR
  return(metaResDf)
}

metaResDf = getResBind(path = "../Vivacilty_v1.0.1/process/", freqCol = 'ALLELE.FREQUENCY', minFreq = 0, maxFreq = 1)
metaResDf = addVarSB(metaResDf)

write.csv(metaResDf, "test_consensus/metaseqLofreq.csv", row.names = F)