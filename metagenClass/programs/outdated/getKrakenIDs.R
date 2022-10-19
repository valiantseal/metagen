

setwd('/home/ubuntu/extraVol/metagenClass/2022-02-25/fastpKrakUniq/kraqSummary')



virusId<-function(x){
  filesList<-list.files(pattern = '.report')
  virusTable<-data.frame(matrix(ncol=0, nrow=1))
  for ( i in filesList){
    sumVirus<-data.frame(matrix(ncol=0, nrow=1))
    df<-read.delim(i, comment.char="#")
    virus<-as.data.frame( df[grep(x,  df$taxName, ignore.case=T),])
    virusID<-data.frame(virus$taxID)
    colnames(virusID)<-'ID'
    virusTable<-rbind(virusTable,virusID)
  }
  fileName<-gsub(" ", "_", x)
  uniqID<-data.frame(unique(virusTable$ID))
  colnames(uniqID)<-'krakenID'
  uniqID$Virus<-x
  #write.table(uniqID, paste0(fileName, ".krakenid"), row.names = F, col.names = F, quote = F)
  return(uniqID)
}

mast<-virusId(x='Human mastadenovirus D')

alpha<-virusId(x='Human alphaherpesvirus 2')

poly<-virusId(x='JC polyomavirus')

combinedVirus<-rbind(mast, alpha, poly)

idList<-data.frame(combinedVirus$krakenID)

write.table(idList, 'kraken.ids', row.names = F, col.names = F, quote = F)

#
krakenReads<-list.files(pattern = '.reads')

#i<-read.delim('B2E22-002A_S1_L001.reads', F, sep = '\t')
#i_filt<-i[(i$V3%in%combinedVirus$krakenID),]
#selCol<-i_filt[, c('V2', 'V3')]
#colnames(selCol)<-c('Read', 'krakenID')
#dfVir<-plyr::join(selCol, combinedVirus, by='krakenID', type='left', match='first')


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
write.csv(allReadsSum, './combined/krakUniqSelVirResultsUniq.csv', row.names = F)

write.table(allKrakenReads, './combined/krakenSelVirReads.tsv', row.names = F, sep = '\t', quote = F)
