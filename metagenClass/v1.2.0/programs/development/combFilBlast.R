dirs<-read.table('newdir.list')

dirs<-as.character(dirs$V1)

combBlast<-function(){
  
  blastComb<-data.frame(matrix(nrow = 0, ncol = 0))
  for (dir in dirs){
    targSample<-paste0('./process/', dir, '/splitSeq10K/')
    filesList<-list.files(targSample)
    for (i in filesList){
      targFile<-paste0(targSample, i, '/blastResFilt.tsv')
      if (file.exists(targFile) == T) {
        curDf<-read.delim(targFile, header = T, sep = '\t')
        curDf$Sample<-dir
        colnames(curDf)[1]<-'Read'
        blastComb<-unique(rbind(blastComb, curDf))
        print(paste0(targFile, '_______done'))
      }
    }
  }
  blastCombUn<-unique(blastComb)
  return(blastCombUn)
}

system.time({blastCombUn<-combBlast()})


write.table(blastCombUn, './blastNtSummary/blastNtSelVirReads.tsv', row.names = F, sep = '\t', quote = F)

resSum<-data.frame(table(blastCombUn$Sample, blastCombUn$Virus))

write.csv(resSum, './blastNtSummary/blastNtSelVirSummary.csv', row.names = F)