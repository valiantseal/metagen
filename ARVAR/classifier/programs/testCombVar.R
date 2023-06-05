library(ggplot2)
library(ggrepel)

testDat = read.csv('test_consensus/metaseqLofreq.csv')

refDat = read.csv("Ludy_metaAmpIvar_overlapSnv.csv")

refFilt = refDat[(refDat$ALT_FREQ >= 0.02) & (refDat$ConsTest == 1),]
testFilt = testDat[(testDat$ALLELE.FREQUENCY >= 0.02),]

length(unique(testFilt$Samp_Pos_Ref_Alt))
length(unique(refFilt$Samp_Pos_Ref_Alt))

length(unique(testFilt$Samp_Pos_Ref_Alt[testFilt$Samp_Pos_Ref_Alt%in%refFilt$Samp_Pos_Ref_Alt]))

allConsSnv = unique(refFilt$Samp_Pos_Ref_Alt)


varAlPos = c(0.2, 0.3, 0.4)
sb = c(11,13,23,37, 53)
# overall summary function
getSummary = function(df, trueSet, varAlPos, sb) {
  allSampleName = character()
  allSnv = numeric()
  allTrue = numeric()
  allTruePos = numeric()
  allFalsePos = numeric()
  allVarPos = numeric()
  allSb = numeric()
  dfSample = df
  for (curVarAlPos in varAlPos) {
    for (curSb in sb) {
      dfFilt = dfSample[(df$Var_Al_RelPos >= curVarAlPos) & (df$STRAND.BIAS <= curSb),]
      curPosNumber = length(unique(dfFilt$Samp_Pos_Ref_Alt[(dfFilt$Samp_Pos_Ref_Alt%in%trueSet)]))
      curFalseNumber = length(unique(dfFilt$Samp_Pos_Ref_Alt[!(dfFilt$Samp_Pos_Ref_Alt%in%trueSet)]))
      curTotNumber = length(unique(dfFilt$Samp_Pos_Ref_Alt))
      curTrueNumb = length(trueSet)
      
      allSnv = c(allSnv, curTotNumber)
      allTrue = c(allTrue, curTrueNumb)
      allTruePos = c(allTruePos, curPosNumber)
      allFalsePos = c(allFalsePos, curFalseNumber)
      allVarPos = c(allVarPos, curVarAlPos)
      allSb = c(allSb, curSb)
    }
  }
  combDat = data.frame( allVarPos, allSb, allSnv, allTrue, allTruePos, allFalsePos)
  colnames(combDat) = c("Var_Al_Pos", "Strand_Bias", "Total_SNVs", "Total_Truth", "True_Positive", "False_Positive")
  return(combDat)
}

allSum = getSummary(df = testFilt, trueSet = allConsSnv, varAlPos = varAlPos, sb = sb)

p2 =  ggplot(allSum, aes(x = False_Positive, y = True_Positive)) +
  geom_text(aes(label=Var_Al_Pos)) + 
  geom_label_repel(aes(label = Strand_Bias),
                   box.padding   = 0.25, 
                   point.padding = 0.25,
                   segment.color = 'grey50', max.overlaps = 50) +
  ggtitle("All samples and libraries combined") +
  theme(plot.title = element_text(hjust = 0.5, vjust = 2))
p2


# per sample summary 
getPerSampSummary = function(df, trueSet, varAlPos, sb) {
  allSampleName = character()
  allSnv = numeric()
  allTrue = numeric()
  allTruePos = numeric()
  allFalsePos = numeric()
  allVarPos = numeric()
  allSb = numeric()
  
  samplesList = unique(df$Sample)
  
  for (sample in samplesList) {
    dfSample = df[(df$Sample == sample),]
    sampleName = unique(dfSample$Sample)
    curTrue = grep(sampleName, allConsSnv, value = T) 
    curTrueNumb = length(curTrue)
    for (curVarAlPos in varAlPos) {
      for (curSb in sb) {
        dfFilt = dfSample[(df$Var_Al_RelPos >= curVarAlPos) & (df$STRAND.BIAS <= curSb),]
        curPosNumber = length(unique(dfFilt$Samp_Pos_Ref_Alt[(dfFilt$Samp_Pos_Ref_Alt%in%curTrue)]))
        curFalseNumber = length(unique(dfFilt$Samp_Pos_Ref_Alt[!(dfFilt$Samp_Pos_Ref_Alt%in%curTrue)]))
        curTotNumber = length(unique(dfFilt$Samp_Pos_Ref_Alt))
        
        allSampleName = c(allSampleName, sample)
        allSnv = c(allSnv, curTotNumber)
        allTrue = c(allTrue, curTrueNumb)
        allTruePos = c(allTruePos, curPosNumber)
        allFalsePos = c(allFalsePos, curFalseNumber)
        allVarPos = c(allVarPos, curVarAlPos)
        allSb = c(allSb, curSb)
      }
    }
  }
  combDat = data.frame(allSampleName, allVarPos, allSb, allSnv, allTrue, allTruePos, allFalsePos)
  colnames(combDat) = c('Sample', "Var_Al_Pos", "Strand_Bias", "Total_SNVs", "Total_Truth", "True_Positive", "False_Positive")
  return(combDat)
}

sumResSamp = getPerSampSummary(df = testFilt, trueSet = allConsSnv, varAlPos = varAlPos, sb = sb)


p3 =  ggplot(sumResSamp, aes(x = False_Positive, y = True_Positive)) +
  geom_text(aes(label=Var_Al_Pos)) + 
  geom_label_repel(aes(label = Strand_Bias),
                   box.padding   = 0.25, 
                   point.padding = 0.25,
                   segment.color = 'grey50', max.overlaps = 50) +
  ggtitle("Per Sample") +
  theme(plot.title = element_text(hjust = 0.5, vjust = 2))

p3

# save results 
targDir = "test_consensus/combVAr/"
dir.create(targDir, showWarnings = F, recursive = T)

write.csv(allSum, paste0(targDir, "LudyMetaseq_Vivacity_VarAlRP_SB_AllComb.csv"), row.names = F)
write.csv(sumResSamp, paste0(targDir, "LudyMetaseq_Vivacity_VarAlRP_SB_PerSample.csv"), row.names = F)

ggsave(filename = paste0(targDir, "LudyMetaseq_Vivacity_VarAlRP_SB_AllComb.png"), plot = p2, height = 10, width = 12, units = "in", dpi = 300)
#ggsave(filename = paste0(targDir, "LudyMetaseq_Vivacity_VarAlRP_SB_PerSample_1.png"), plot = p3, height = 10, width = 12, units = "in", dpi = 300)