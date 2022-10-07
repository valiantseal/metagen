

setwd('/home/ubuntu/extraVol/metagenClass/2022-02-25/fastpKrakUniq/kraqSummary')



selVirus<-function(x){
  filesList<-list.files(pattern = '.report')
  virusTable<-data.frame(matrix(ncol=0, nrow=1))
  for ( i in filesList){
    sumVirus<-data.frame(matrix(ncol=0, nrow=1))
    df<-read.delim(i, comment.char="#")
    virus<-as.data.frame( df[grep(x,  df$taxName, ignore.case=T),])
    sumVirus$Sample<-i
    sumVirus$Virus<-x
    if (nrow(virus) > 0 ){
      sumVirus$Reads<-sum(virus$reads)
    } else {
      sumVirus$Reads<-0
    }
    #combReads<-data.frame(matrix(ncol=0, nrow=1))
    #combReads$
    
    
    virusTable<-rbind(virusTable, sumVirus)
    virusTable$Sample<-gsub('.report', '', virusTable$Sample)
  }
  return(virusTable)
}

mastedSum<-selVirus(x='Human mastadenovirus D')

alphaSum<-selVirus(x='Human alphaherpesvirus 2')

polySum<-selVirus(x='JC polyomavirus')

combSum<-rbind(alphaSum, mastedSum, polySum)

write.csv(combSum, './combined/krakUniqSelVirResults.csv', row.names = F)

# kraken extract read

