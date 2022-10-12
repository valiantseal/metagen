setwd('/home/ubuntu/extraVol/metagenClass/2022-02-25/fastpKrakUniq/kraqSummary')

combinedVirus<-read.csv('summaryKrakenIDs.csv')
krakenReads<-list.files(pattern = '.reads')

allKrakenReads<-data.frame(matrix(ncol = 0, nrow = 0))

for (sample in krakenReads){
  i<-read.delim(sample, F, sep = '\t')
  i_filt<-i[(i$V3%in%combinedVirus$krakenID),]
  selCol<-i_filt[, c('V2', 'V3')]
  colnames(selCol)<-c('Read', 'krakenID')
  dfVir<-plyr::join(selCol, combinedVirus, by='krakenID', type='left', match='first')
  dfVirUniq<-unique(dfVir[, c('Read', 'Virus')])
  dfVirUniq$Sample<-sample
  allKrakenReads<-rbind(allKrakenReads, dfVirUniq)
}

allKrakenReads$Sample<-gsub('.reads', '', allKrakenReads$Sample)

allReadsSum<-data.frame(table(allKrakenReads$Sample, allKrakenReads$Virus))
colnames(allReadsSum)[1:2]<-c('Sample', 'Virus')


write.csv(allReadsSum, 'krakTargetVirSummary.csv', row.names = F)


# check that all reads identyfy only with one virus
idCheck<-data.frame(table(allKrakenReads$Sample, allKrakenReads$Read))

extrCh<-allKrakenReads[(allKrakenReads$Read=='M01541:32:000000000-K7WN6:1:1101:12606:18147'),]

write.table(allKrakenReads, 'krakenSelVirReads.tsv', row.names = F, sep = '\t', quote = F)




