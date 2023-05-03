samplesList = c("2884X", "2885Y", "2886Z")


# get results per sample
getRes = function(path, pat) {
  combTab = list()
  filesList = list.files(path, pattern = pat)
  for ( i in filesList) {
    inFile =  paste0(path, i)
    df = read.delim(inFile, T)
    dfFilt = df[(df$PASS== "TRUE") & (df$ALT_FREQ >= 0.02),]
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

ampSnvs = unique(getConsAllSamp(samplesList = samplesList, path = "test_consensus/viralrecon_ampseq/"))
length(ampSnvs)

metaSnvs = unique(getConsAllSamp(samplesList = samplesList, path = "test_consensus/viralrecon_metaseq/"))
length(unique(metaSnvs))

allConsSnv = unique(metaSnvs[(metaSnvs%in%ampSnvs)])
length(allConsSnv)

# get vivacity data
metaResDf = read.csv("test_consensus/vivacity_metaConsCheck.csv")

metaResDf$ConsTest = NA
metaResDf$ConsTest[!(metaResDf$Samp_Pos_Ref_Alt %in% allConsSnv)] = 0
metaResDf$ConsTest[(metaResDf$Samp_Pos_Ref_Alt %in% allConsSnv)] = 1

sumConsTest = data.frame(table(metaResDf$ConsTest))

colnames(metaResDf)

# models

# multivariate
multivarModel = glm(ConsTest ~ ALLELE.FREQUENCY  + STRAND.BIAS + DEPTH + QUAL + Ref_Al_RelPos + Var_Al_Relpos, data = metaResDf, family = binomial)
summary(multivarModel)

multivarModel = glm(ConsTest ~ ALLELE.FREQUENCY  + STRAND.BIAS + Ref_Al_RelPos + Var_Al_Relpos, data = metaResDf, family = binomial)
summary(multivarModel)


# multivariate
multivarModel = glm(ConsTest ~ ALLELE.FREQUENCY  + STRAND.BIAS + DEPTH + QUAL + Position_test, data = metaResDf, family = binomial)
summary(multivarModel)

# run all univar models
vars = c('ALLELE.FREQUENCY', 'STRAND.BIAS',  'DEPTH', 'Ref_Al_RelPos', "QUAL", 'Var_Al_Relpos', "Position_test", "Freq_adj", "Type", "Pi.Ln.Pi.")

runUnivarMod = function(vars, df) {
  resultsDf  = data.frame(matrix(nrow = 0, ncol = 0))
  for (curVar in vars) {
    testVar = df[["ConsTest"]]
    indVar = df[[curVar]]
    univarModel = glm(testVar ~ indVar, family = binomial)
    sumUniv = summary(univarModel)
    sumUnivDf = data.frame(sumUniv$coefficients)
    dfSub =  sumUnivDf[2,]
    dfSub$Var = curVar
    resultsDf = rbind(resultsDf, dfSub)
  }
  return(resultsDf)
}

univarRes = runUnivarMod(vars = vars, df = metaResDf)
colnames(univarRes)[4] = "Pval"

#write.csv(univarRes, file = "test_consensus/vivacity_ivar_metaConsCheck_univarRes.csv", row.names = F)


mean(metaResDf$STRAND.BIAS[metaResDf$ConsTest == 1])
median(metaResDf$STRAND.BIAS[metaResDf$ConsTest == 1])
