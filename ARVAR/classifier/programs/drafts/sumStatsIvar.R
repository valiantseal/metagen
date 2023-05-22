samplesList = c("2884X", "2885Y", "2886Z")

# get vivacity data
getResBind = function(path, pat, freqCol, minFreq, maxFreq) {
  combTab = data.frame(matrix(nrow = 0 , ncol = 0))
  filesList = list.files(path, pattern = pat)
  for ( i in filesList) {
    inFile = paste0(path, i, "/filtered.csv")
    df = read.csv(inFile)
    df$RawVarFreq = (df$FWD.VAR + df$REV.VAR) / (df$FWD.VAR + df$REV.VAR + df$FWD.REF + df$REV.REF)
    dfFilt = df[(df[freqCol] >= minFreq) & (df[freqCol] <= maxFreq),]
    dfFilt$Sample = pat
    dfFilt$Samp_Pos_Ref_Alt = paste(dfFilt$Sample, dfFilt$POSITION, dfFilt$REF.NT, dfFilt$VAR.NT, sep = "__")
    
    combTab = rbind(combTab, dfFilt)
  }
  return(combTab)
}

# get vivacity results for all samples
getAllResBind = function(samplesList, path, freqCol, minFreq, maxFreq) {
  allSnv = data.frame(matrix(nrow = 0 , ncol = 0))
  for (sample in samplesList) {
    dfList = getResBind(path = path, pat = sample, freqCol = freqCol, minFreq = minFreq, maxFreq = maxFreq)
    allSnv = rbind(allSnv, dfList)
  }
  return(allSnv)
}

# calculate strand bias ratio for variant 
addVarSB = function(metaResDf) {
  metaResDf$Var_SB = metaResDf$FWD.VAR / metaResDf$REV.VAR
  return(metaResDf)
}

metaResDf = getAllResBind(samplesList = samplesList, path = "process/", freqCol = 'ALLELE.FREQUENCY', minFreq = 0, maxFreq = 1)
metaResDf = addVarSB(metaResDf)

## get ivar results
# get results per sample
getRes = function(path, pat) {
  combTab = list()
  filesList = list.files(path, pattern = pat)
  for ( i in filesList) {
    inFile =  paste0(path, i)
    df = read.delim(inFile, T)
    dfFilt = df[(df$PASS== "TRUE") & (df$ALT_FREQ >= 0.01),]
    dfFilt$Sample = pat
    dfFilt$Samp_Pos_Ref_Alt = paste(dfFilt$Sample, dfFilt$POS, dfFilt$REF, dfFilt$ALT, sep = "__")
    
    combTab = append(combTab, list(dfFilt))
  }
  return(combTab)
}


getConsensus = function(dfList) {
  df1 = dfList[[1]]$Samp_Pos_Ref_Alt
  df2 = dfList[[2]]$Samp_Pos_Ref_Alt
  combDf = df1[(df1%in%df2)]
  return(combDf)
} 

getConsAllSamp = function(samplesList, path) {
  allSnv = character()
  for (sample in samplesList) {
    dfList = getRes(path = path, pat = sample)
    curSnv = getConsensus(dfList =  dfList)
    allSnv = c(allSnv, curSnv)
  }
  return(allSnv)
}

ampSnvs = unique(getConsAllSamp(samplesList = samplesList, path = "../test_consensus/viralrecon_ampseq/"))
length(ampSnvs)

metaSnvs = unique(getConsAllSamp(samplesList = samplesList, path = "../test_consensus/viralrecon_metaseq/"))
length(unique(metaSnvs))

allConsSnv = unique(metaSnvs[(metaSnvs%in%ampSnvs)])
length(allConsSnv)

metaResDf$ConsTest = NA
metaResDf$ConsTest[!(metaResDf$Samp_Pos_Ref_Alt %in% allConsSnv)] = 0
metaResDf$ConsTest[(metaResDf$Samp_Pos_Ref_Alt %in% allConsSnv)] = 1

sumConsTest = data.frame(table(metaResDf$ConsTest))

length(unique(metaResDf$Samp_Pos_Ref_Alt))
length(unique(metaResDf$Samp_Pos_Ref_Alt[metaResDf$Samp_Pos_Ref_Alt%in%allConsSnv]))
length(allConsSnv[allConsSnv%in%metaResDf$Samp_Pos_Ref_Alt])

write.csv(metaResDf, "test_consensus/vivacityNoFilt_ivarfreq0.01.csv", row.names = F)