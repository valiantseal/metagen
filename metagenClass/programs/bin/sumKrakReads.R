setwd('./kraqSummary')

# read kraken IDs present in our samples
combinedVirus<-read.csv('summaryKrakenIDs.csv')

# list files with reads from kraken
krakenReads<-list.files(pattern = '.reads')

allKrakenReads<-data.frame(matrix(ncol = 0, nrow = 0))

# for each sample add the name of the virus to the read from kraken, make sure reads/virus combinations are unique
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

# edit sample names
allKrakenReads$Sample<-gsub('.reads', '', allKrakenReads$Sample)

# summarize number of reads for each sample/virus combination
allReadsSum<-data.frame(table(allKrakenReads$Sample, allKrakenReads$Virus))
colnames(allReadsSum)[1:2]<-c('Sample', 'Virus')

write.csv(allReadsSum, 'krakTargetVirSummary.csv', row.names = F)

write.table(allKrakenReads, 'krakenSelVirReads.tsv', row.names = F, sep = '\t', quote = F)


# check that all reads identyfy only with one virus
idCheck<-data.frame(table(allKrakenReads$Sample, allKrakenReads$Read))

extrCh<-allKrakenReads[(allKrakenReads$Read=='M01541:32:000000000-K7WN6:1:1101:12606:18147'),]





