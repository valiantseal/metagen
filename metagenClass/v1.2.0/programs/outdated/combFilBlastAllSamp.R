dirs<-read.table('newdir.list')

dirs<-as.character(dirs$V1)

combBlast<-function(){
  
  blastComb<-data.frame(matrix(nrow = 0, ncol = 0))
  for (dir in dirs){
    targSample<-paste0('./process/', dir, '/blastFiltComb.tsv')
    if (file.exists(targSample) == T) {
      curDf<-read.delim(targSample, header = T, sep = '\t')
      blastComb<-rbind(blastComb, curDf)
      print(paste0(targSample, '_______done'))
    }
  }
  blastCombUn<-unique(blastComb)
  return(blastCombUn)
}




blastCombUn<-combBlast()


write.table(blastCombUn, './blastNtSummary/blastNtSelVirReads.tsv', row.names = F, sep = '\t', quote = F)

resSum<-data.frame(table(blastCombUn$Sample, blastCombUn$Virus))

write.csv(resSum, './blastNtSummary/blastNtSelVirSummary.csv', row.names = F)