setwd('./kraqSummary')



virusId<-function(x){
  filesList<-list.files(pattern = '.report')
  virusTable<-data.frame(matrix(ncol=0, nrow=1))
  # for each kraken sample file extract rows with matching viruses
  for ( i in filesList){
    sumVirus<-data.frame(matrix(ncol=0, nrow=1))
    df<-read.delim(i, comment.char="#")
    virus<-as.data.frame( df[grep(x,  df$taxName, ignore.case=T),])
    virusID<-data.frame(virus$taxID)
    colnames(virusID)<-'ID'
    virusTable<-rbind(virusTable,virusID)
  }
  # make sure function works if no tarket virus is found and keep only unique kraken IDs
  if (nrow(virusTable)>0){
    uniqID<-data.frame(unique(virusTable$ID))
    colnames(uniqID)<-'krakenID'
  } else {
    uniqID<-data.frame('No_reads')
    colnames(uniqID)<-'krakenID'
  }
  uniqID$Virus<-x
  return(uniqID)
}

# read list of target viruses
virusList<-read.table('../virus.list', F, sep='\t')
virus.list<-virusList$V1

# bind kraken id's for all viruses from the list
sumStat<-data.frame(matrix(ncol=0, nrow=0))
for (virus in virus.list){
  pattern<-gsub(" ", "_", virus)
  sumVirus<-virusId(x=virus)
  sumStat<-rbind(sumStat, sumVirus)
}

combinedVirus<-sumStat[!(sumStat$krakenID=='No_reads'),]

idList<-data.frame(combinedVirus$krakenID)

write.table(idList, 'kraken.ids', row.names = F, col.names = F, quote = F)

write.csv(sumStat, 'summaryKrakenIDs.csv', row.names = F)