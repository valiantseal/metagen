combResult = function(path) {
  combData = data.frame(matrix(nrow = 0, ncol = 0))
  filesList = list.files(path)
  for ( i in filesList) {
    try({
      df = read.delim(paste0(path,i))
      dfFilter = df[(df$PASS == "TRUE"),]
      sample = gsub("\\..*", "", i)
      dfFilter$FullSamp = sample
      sample = gsub("\\_.*", "", sample)
      dfFilter$ExactSamp = sample
      sampleList = strsplit(sample, '-')
      sampleName =  sampleList[[1]][1:3]
      sampleName = paste(sampleName, collapse = "-")
      dfFilter$Sample = sampleName
      dfFilter$Samp_Pos_Ref_Alt = paste(dfFilter$Sample, dfFilter$POS, dfFilter$REF, dfFilter$ALT, sep = "__")
      combData = rbind(combData, dfFilter)
    })
  }
  return(combData )
}


getRepeating = function(testList) {
  unList = unique(testList)
  newList = character()
  for (i in unList) {
    sample = gsub("\\_.*", "", i)
    sampleList = strsplit(sample, '-')
    sampleName =  sampleList[[1]][1:3]
    sampleName = paste(sampleName, collapse = "-")
    newList = c(newList, sampleName)
  }
  sumDat = data.frame(table(newList))
  colnames(sumDat)[1] = "Sample"
  return(sumDat)
}

metaSeq = combResult(path = "/home/ubuntu/extraVol/Viralrecon/covid/Ludy_metaseq_2023-05-16/output/relaxIvar/")
length(unique(metaSeq$Sample))
ampSeq = combResult(path = "/home/ubuntu/extraVol/Viralrecon/covid/Ludy_Apr242023/output/relaxIvar/")
length(unique(ampSeq$Sample))

length(unique(metaSeq$Sample[metaSeq$Sample%in%ampSeq$Sample]))


# find local consensus 
ampSum = getRepeating(testList = ampSeq$FullSamp)
colnames(ampSum)

getLocCons = function(sumDat, testDat) {
  singleList = as.character(unique(sumDat$Sample[(sumDat$Freq < 2)]))
  doubList =  as.character(unique(sumDat$Sample[(sumDat$Freq > 1)]))
  doubCons = character()
  for (i in doubList) {
    subDf = testDat[(testDat$Sample == i),]
    sumSab = data.frame(table(subDf$Samp_Pos_Ref_Alt))
    colnames(sumSab)[1] = "Samp_Pos_Ref_Alt"
    conSnv = as.character(unique(sumSab$Samp_Pos_Ref_Alt[(sumSab$Freq > 1)]))
    doubCons = c(doubCons, conSnv)
  }
  doubConsUniq = unique(doubCons)
  singDat = testDat[(testDat$Sample %in% singleList),]
  doubDat = testDat[(testDat$Sample %in% doubList),]
  doubFilt = doubDat[(doubDat$Samp_Pos_Ref_Alt %in% doubConsUniq),]
  combDat = rbind(singDat, doubFilt)
  return(combDat)
}

ampSeqCons = getLocCons(sumDat = ampSum, testDat = ampSeq)
length(unique(ampSeq$Sample))
length(unique(ampSeqCons$Sample))

dir.create("test_consensus")

write.csv(metaSeq, "test_consensus/metaseqIvar.csv", row.names = F)
write.csv(ampSeq, "test_consensus/ampseqIvar.csv", row.names = F)
write.csv(ampSeqCons, "test_consensus/ampseqConsIvar.csv", row.names = F)