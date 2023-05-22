library(ggrepel)
library(psych)

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

metaResDf = getAllResBind(samplesList = samplesList, path = "process/", freqCol = 'RawVarFreq', minFreq = 0.02, maxFreq = 0.98)
metaResDf = addVarSB(metaResDf)

metaResDf$ConsTest = NA
metaResDf$ConsTest[!(metaResDf$Samp_Pos_Ref_Alt %in% allConsSnv)] = 0
metaResDf$ConsTest[(metaResDf$Samp_Pos_Ref_Alt %in% allConsSnv)] = 1

sumConsTest = data.frame(table(metaResDf$ConsTest))

colnames(metaResDf)


# get summary numbers 

# Overlap between true set and current set
length(unique(metaResDf$Samp_Pos_Ref_Alt))
length(unique(metaResDf$Samp_Pos_Ref_Alt[metaResDf$Samp_Pos_Ref_Alt%in%allConsSnv]))
# SNPs in Daras but not mine
allConsSnv[!(allConsSnv%in%metaResDf$Samp_Pos_Ref_Alt)]

# work with filtered Dara's snvs with particular raw frequencies
source('~/extraVol/ARVAR/Vivacity/Vivacilty_v1.0.1/programs/drafts/getPrevLofreqRes.R')
length(unique(metaResDf$Samp_Pos_Ref_Alt[metaResDf$Samp_Pos_Ref_Alt%in%selConsSnv]))
selConsSnv[!(selConsSnv%in%metaResDf$Samp_Pos_Ref_Alt)]

metaResFilt = metaResDf[(metaResDf$Type == "Substitution"),]
length(unique(metaResFilt$Samp_Pos_Ref_Alt))
length(unique(metaResFilt$Samp_Pos_Ref_Alt[metaResFilt$Samp_Pos_Ref_Alt%in%selConsSnv]))

##
# test different variables


findOptim = function(df, test_var, test_val, allConsSnv, smaller) {
  allSnps = numeric()
  allTPos = numeric()
  allFPos = numeric()
  allRecall = numeric()
  allFPR = numeric()
  allVals = numeric()
  allPsedAUC = numeric()
  allHarmMean = numeric()
  Ratio = numeric()
  for (curVal in test_val) {
    if (smaller == T) {
      dfFilt = df[(df[[test_var]] <= curVal),]
    } else {
      dfFilt = df[(df[[test_var]] >= curVal),]
    }
    curSnps = length(unique(dfFilt$Samp_Pos_Ref_Alt))
    allSnps = c(allSnps, curSnps)
    curTPos = length(unique(dfFilt$Samp_Pos_Ref_Alt[(dfFilt$Samp_Pos_Ref_Al %in% allConsSnv)]))
    allTPos = c(allTPos, curTPos)
    curFPos = length(unique(dfFilt$Samp_Pos_Ref_Alt[!(dfFilt$Samp_Pos_Ref_Al %in% allConsSnv)]))
    allFPos = c(allFPos, curFPos)
    curRecall = curTPos / length(allConsSnv)
    allRecall = c(allRecall, curRecall)
    curFPR = curFPos / curSnps
    allFPR = c(allFPR, curFPR)
    allVals = c(allVals, curVal)
    curHarmMean = harmonic.mean(curRecall, curFPR)
    curAUC = 1/2 - curRecall/2 + curFPR/2
    allHarmMean = c(allHarmMean, curHarmMean)
    allPsedAUC = c(allPsedAUC, curAUC)
    curRatio = curRecall / curFPR
    Ratio = c(Ratio, curRatio)
  }
  combData = data.frame(allSnps, allTPos, allFPos, allRecall, allFPR, allHarmMean, allPsedAUC, Ratio,  allVals)
  colnames(combData) = c("SNPs", "Tpos", "Fpos", "Recall", "FPR", "HarmonicMean", "PseudoAUC", "Ratio", 'VarNumb')
  return(combData)
}


makePlot = function(df, test_var, test_val, allConsSnv, smaller) {
  datFilt = df[!(df[test_var] == Inf),]
  combDat = findOptim(df = datFilt , test_var = test_var, test_val = test_val, allConsSnv = allConsSnv, smaller = smaller)
 p1 =  ggplot(combDat, aes(x = Fpos, y = Tpos)) +
    geom_point(size = 0.1) + 
    geom_label_repel(aes(label = VarNumb),
                     box.padding   = 0.25, 
                     point.padding = 0.25,
                     segment.color = 'grey50', max.overlaps = 50) +
    ggtitle(test_var) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 2))
 print(p1)
 ggsave(filename = paste0("test_consensus/freqFil_", test_var, ".png"), plot = p1, width = 14, height = 10, units = "in", dpi = 300)
}

colnames(metaResFilt)
# strandBias
test_val = seq(from = 1, by = 2, length.out = 100) 
makePlot(df = metaResFilt, test_var = "STRAND.BIAS", test_val = test_val, allConsSnv = selConsSnv, smaller = T)
test_val = seq(from = 0, by = 0.2, length.out = 100) 
makePlot(df = metaResFilt, test_var = "Var_SB", test_val = test_val, allConsSnv = selConsSnv, smaller = T)

# variant relative position
test_val = seq(from = 0, by = 0.01, length.out = 100) 
makePlot(df = metaResFilt, test_var = "Var_Al_RelPos", test_val = test_val, allConsSnv = selConsSnv, smaller = F)


## try comparisons without filtering Dara's dataset
metaResDf = getAllResBind(samplesList = samplesList, path = "process/", freqCol = 'RawVarFreq', minFreq = 0.02, maxFreq = 1)
metaResDf = addVarSB(metaResDf)
metaResFilt = metaResDf[(metaResDf$Type == "Substitution"),]
length(unique(metaResFilt$Samp_Pos_Ref_Alt))
length(unique(metaResFilt$Samp_Pos_Ref_Alt[metaResFilt$Samp_Pos_Ref_Alt%in%allConsSnv]))

makePlot = function(df, test_var, test_val, allConsSnv, smaller) {
  datFilt = df[!(df[test_var] == Inf),]
  combDat = findOptim(df = datFilt , test_var = test_var, test_val = test_val, allConsSnv = allConsSnv, smaller = smaller)
  p1 =  ggplot(combDat, aes(x = Fpos, y = Tpos)) +
    geom_point(size = 0.1) + 
    geom_label_repel(aes(label = VarNumb),
                     box.padding   = 0.25, 
                     point.padding = 0.25,
                     segment.color = 'grey50', max.overlaps = 50) +
    ggtitle(test_var) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 2))
  print(p1)
  ggsave(filename = paste0("test_consensus/freqFilMinOnly_", test_var, ".png"), plot = p1, width = 14, height = 10, units = "in", dpi = 300)
}

# strandBias
test_val = seq(from = 1, by = 2, length.out = 100) 
makePlot(df = metaResFilt, test_var = "STRAND.BIAS", test_val = test_val, allConsSnv = allConsSnv, smaller = T)
test_val = seq(from = 0, by = 0.2, length.out = 100) 
makePlot(df = metaResFilt, test_var = "Var_SB", test_val = test_val, allConsSnv = allConsSnv, smaller = T)

# variant relative position
test_val = seq(from = 0, by = 0.01, length.out = 100) 
makePlot(df = metaResFilt, test_var = "Var_Al_RelPos", test_val = test_val, allConsSnv = allConsSnv, smaller = F)
makePlot(df = metaResFilt, test_var = "Ref_Al_RelPos", test_val = test_val, allConsSnv = allConsSnv, smaller = F)


# no filtering  compare frequencies
## try comparisons without filtering Dara's dataset
metaResDf = getAllResBind(samplesList = samplesList, path = "process/", freqCol = 'RawVarFreq', minFreq = 0, maxFreq = 1)
metaResDf = addVarSB(metaResDf)
metaResFilt = metaResDf[(metaResDf$Type == "Substitution"),]
length(unique(metaResFilt$Samp_Pos_Ref_Alt))
length(unique(metaResFilt$Samp_Pos_Ref_Alt[metaResFilt$Samp_Pos_Ref_Alt%in%allConsSnv]))

makePlot = function(df, test_var, test_val, allConsSnv, smaller) {
  datFilt = df[!(df[test_var] == Inf),]
  combDat = findOptim(df = datFilt , test_var = test_var, test_val = test_val, allConsSnv = allConsSnv, smaller = smaller)
  p1 =  ggplot(combDat, aes(x = Fpos, y = Tpos)) +
    geom_point(size = 0.1) + 
    geom_label_repel(aes(label = VarNumb),
                     box.padding   = 0.25, 
                     point.padding = 0.25,
                     segment.color = 'grey50', max.overlaps = 50) +
    ggtitle(test_var) +
    theme(plot.title = element_text(hjust = 0.5, vjust = 2))
  print(p1)
  ggsave(filename = paste0("test_consensus/NoFiltering_", test_var, ".png"), plot = p1, width = 14, height = 10, units = "in", dpi = 300)
}

colnames(metaResFilt)
# strandBias
test_val = seq(from = 1, by = 2, length.out = 100) 
makePlot(df = metaResFilt, test_var = "STRAND.BIAS", test_val = test_val, allConsSnv = allConsSnv, smaller = T)
test_val = seq(from = 0, by = 0.2, length.out = 100) 
makePlot(df = metaResFilt, test_var = "Var_SB", test_val = test_val, allConsSnv = allConsSnv, smaller = T)

# variant relative position
test_val = seq(from = 0, by = 0.01, length.out = 100) 
makePlot(df = metaResFilt, test_var = "Var_Al_RelPos", test_val = test_val, allConsSnv = allConsSnv, smaller = F)
makePlot(df = metaResFilt, test_var = "Ref_Al_RelPos", test_val = test_val, allConsSnv = allConsSnv, smaller = F)

# frequencies
# variant relative position
test_val = seq(from = 0, by = 0.01, length.out = 100) 
makePlot(df = metaResFilt, test_var = "ALLELE.FREQUENCY", test_val = test_val, allConsSnv = allConsSnv, smaller = F)
makePlot(df = metaResFilt, test_var = "RawVarFreq", test_val = test_val, allConsSnv = allConsSnv, smaller = F)
