df<-read.table('./kraqSummary/krakenSelVirReads.tsv', T, sep = '\t')

samples<-unique(df$Sample)

for ( i in samples){
  sampleDf<-df[(df$Sample==i),]
  readDf<-data.frame(unique(sampleDf$Read))
  filePath=paste0('process/',i,'/krakenReads.id')
  write.table(readDf, file = filePath, row.names = F, col.names = F, quote = F)
}

