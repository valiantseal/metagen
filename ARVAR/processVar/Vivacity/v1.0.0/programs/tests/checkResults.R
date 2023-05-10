filesList = list.files("process")


medStrand = numeric()
meanStran = numeric()
for ( i in filesList) {
  inFile = paste0("process/", i, "/filtered.csv")
  try({
    print(inFile)
    df = read.csv(inFile)
    dfFilt = df[(df$Position_test == "Pass"),]
    curMed = median(dfFilt$STRAND.BIAS)
    curMean = mean(dfFilt$STRAND.BIAS)
    medStrand = c(medStrand, curMed)
    meanStran = c(meanStran , curMean)
  })
}

length(medStrand )
length(meanStran)

median(medStrand)
mean(meanStran)