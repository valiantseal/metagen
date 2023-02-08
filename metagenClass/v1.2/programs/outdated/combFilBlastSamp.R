

combBlast<-function(){
  
  blastComb<-data.frame(matrix(nrow = 0, ncol = 0))
  sampleName<-read.table('sample.name')
  sampleName<-as.character(sampleName$V1)
  targSample<-paste0('./splitSeq10K/')
  filesList<-list.files(targSample)
  for (i in filesList){
    targFile<-paste0(targSample, i, '/blastResFilt.tsv')
    if (file.exists(targFile) == T) {
      curDf<-read.delim(targFile, header = T, sep = '\t')
      curDf$Sample<-sampleName
      colnames(curDf)[1]<-'Read'
      blastComb<-rbind(blastComb, curDf)
      #print(i)
     

    }
  }
  blastCombUn<-unique(blastComb)
  return(blastCombUn)
}

blastCombUn<-combBlast()


write.table(blastCombUn, 'blastFiltComb.tsv', row.names = F, sep = '\t', quote = F)

print(paste0(getwd(), '_______done'))