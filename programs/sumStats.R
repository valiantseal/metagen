samplesList = c("2884X", "2885Y", "2886Z")

refDf = read.csv("test_consensus/iSNVs_final_Dara.csv")

df = refDf
colnames(refDf)
colName = "Allele.frequency.at.Day.0"
sampleName = "2884X"

getTrueSet = function(df, colName, sampleName) {
  refDfFilt = df[!is.na(df[colName]),]
  #refDfFilt = refDfFilt[refDfFilt[, colName] < 0.90,]
  refDfFilt$Samples = sampleName
  refDfFilt$Samp_Pos_Ref_Alt = paste(refDfFilt$Samples, refDfFilt$NT.POSITION, refDfFilt$REF.NT, refDfFilt$VAR.NT , sep = "__")
  snvs = unique(refDfFilt$Samp_Pos_Ref_Alt)
  return(snvs)
}

x = getTrueSet(df = refDf, colName =  "Allele.frequency.at.Day.0" , sampleName = "2884X")
y = getTrueSet(df = refDf, colName =  "Allele.Frequency.at.Day.31" , sampleName = "2885Y")
z = getTrueSet(df = refDf, colName =  "Allele.Frequency.at.Day.44" , sampleName = "2886Z")

allConsSnv = unique(c(x,y,z))
length(allConsSnv)

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

getAllResBind = function(samplesList, path, freqCol, minFreq, maxFreq) {
  allSnv = data.frame(matrix(nrow = 0 , ncol = 0))
  for (sample in samplesList) {
    dfList = getResBind(path = path, pat = sample, freqCol = freqCol, minFreq = minFreq, maxFreq = maxFreq)
    allSnv = rbind(allSnv, dfList)
  }
  return(allSnv)
}

metaResDf = getAllResBind(samplesList = samplesList, path = "process/", freqCol = 'RawVarFreq', minFreq = 0, maxFreq = 1)

metaResDf$ConsTest = NA
metaResDf$ConsTest[!(metaResDf$Samp_Pos_Ref_Alt %in% allConsSnv)] = 0
metaResDf$ConsTest[(metaResDf$Samp_Pos_Ref_Alt %in% allConsSnv)] = 1

sumConsTest = data.frame(table(metaResDf$ConsTest))

colnames(metaResDf)


# get summary numbers 

# Overlap between true set and current set
length(unique(metaResDf$Samp_Pos_Ref_Alt[metaResDf$Samp_Pos_Ref_Alt%in%allConsSnv]))
# SNPs in Daras but not mine
allConsSnv[!(allConsSnv%in%metaResDf$Samp_Pos_Ref_Alt)]