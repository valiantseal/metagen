
combDat = function(path) {
  combDat = data.frame(matrix(nrow = 0, ncol =0))
  samples = list.files(path)
  for (curSample in samples) {
    try({
      curDf = paste0(path, "/", curSample, "/filtered.csv")
      df = read.csv(curDf)
      exactSample = gsub("_", "-", curSample)
      df$ExactSample = exactSample
      sampNamelist = strsplit(exactSample, "-")
      sampleName =  sampNamelist[[1]][1:3]
      sampleName = paste(sampleName, collapse = "-")
      df$Sample = sampleName
      df$Samp_Pos_Ref_Alt = paste(df$Sample, df$POSITION, df$REF.NT, df$VAR.NT, sep = "__")
      
      combDat = rbind(combDat, df)
    })
    
  }
  return(combDat)
}

metaseq = combDat(path = "/home/ubuntu/extraVol/ARVAR/classifier/Vivacilty_v1.0.1/process_par")
ampseq = combDat(path = "/home/ubuntu/extraVol/ARVAR/classifier/Vivacilty_v1.0.1/amp_process_par")

metaseqSamp = unique(metaseq$Sample)
ampFilter = ampseq[(ampseq$Sample %in% metaseqSamp),]

length(unique(metaseq$Sample))
length(unique(ampseq$Sample))
length(unique(ampFilter$Sample))

length(unique(ampFilter$Sample[ampFilter$Sample %in% metaseqSamp]))
length(unique(metaseq$Sample[metaseq$Sample %in% ampFilter$Sample]))

#write.csv(metaseq, "test_consensus/metaseq_comb_vivacity.csv", row.names = F)
#write.csv(ampseq, "test_consensus/ampseq_comb_vivacity.csv", row.names = F)
#write.csv(ampFilter, "test_consensus/ampseq_comb_metaSampFilt_vivacity.csv", row.names = F)

colnames(metaseq)

# find local consensus 


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

ampSum = getRepeating(testList = ampseq$ExactSample)
colnames(ampSum)
