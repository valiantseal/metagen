# read list of target viruses
virusList<-read.table('../../virus.list', F, sep='\t')
virus.list<-virusList$V1

sampleName<-read.table('sample.name')
sampleName<-as.character(sampleName$V1)


virusId<-function(){
  # for each kraken sample file extract rows with matching viruses
  df<-read.delim('krakUniq_sample.report', comment.char="#")
  virus<-as.data.frame(df[grep(paste(virus.list,collapse="|"), df$taxName, ignore.case=T),])
  virus$Virus<-NA
  if (nrow(virus) > 0){
    for ( i in virus.list){
      virus$Virus[grepl(i, virus$taxName, ignore.case = T)]<-i
      }
  }
  virusID<-data.frame(unique(virus[, c('taxID', 'Virus')]))
  colnames(virusID)[1]<-'krakId'
  return(virusID)
}


idList<-virusId()

#write.table(idList, 'kraken.ids', row.names = F, col.names = F, quote = F)



readNameTab<-read.delim('krakUniq_sample.kraken', F, sep = '\t')
colnames(readNameTab)[3]<-'krakId'
colnames(readNameTab)[2]<-'Read'
selCol<-readNameTab[, c('krakId', 'Read')]
# join 

cReads<-plyr::join(idList, selCol, by='krakId', type='left', match='all')

allKrakenReads<-unique(cReads[, c('Virus', 'Read')])
allKrakenReads$Sample<-sampleName

# save sorted kraken output
write.table(allKrakenReads, 'krakenSelVirReads.tsv', row.names = F, sep = '\t', quote = F)

reads<-unique(allKrakenReads[, 'Read'])

write.table(reads, 'krakenReads.id', row.names = F, col.names = F, quote = F)
