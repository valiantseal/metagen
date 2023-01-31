samples<-read.table('newdir.list', F)
samples<-samples$V1

mergeTables<-function(sample) {
  blastIn<-paste0('./process/', sample, '/blastResTop_1.tsv')
  krakIn<-paste0('./process/', sample, '/krakenSelVirReads.tsv')
  
  blastDf<-read.delim(blastIn, T)
  krakDf<-read.delim(krakIn, T)
  blastDf$Sample_Read_Virus<-paste(blastDf$Sample, blastDf$Read, blastDf$Virus, sep = '__')
  krakDf$Sample_Read_Virus<-paste(krakDf$Sample, krakDf$Read, krakDf$Virus, sep = '__')
  
  krakFilt<-krakDf[(krakDf$Sample_Read_Virus%in%blastDf$Sample_Read_Virus),]
  
  filtRes<-unique(krakFilt[, c('Sample', 'Virus', 'Read')])
  
  return(filtRes)
  
}


confirmedRes<-data.frame(matrix(nrow=0, ncol = 0))

for (sample in samples){
  mergeSamp<-mergeTables(sample = sample)
  confirmedRes<-rbind(confirmedRes, mergeSamp)
}

runSummary<-data.frame(table(confirmedRes$Sample, confirmedRes$Virus))

colnames(runSummary)[1:2]<-c('Sample', 'Virus')

write.csv(confirmedRes, 'krakBlastConfReads.csv', row.names = F)
write.csv(runSummary, 'krakBlastConfReads_summary.csv', row.names = F)

