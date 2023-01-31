sumDat<-read.table('./testReads/confirmedReads.list', T)

samples<-unique(sumDat$Sample.x)

filterLength<- 100


for ( i in samples){
  readsPath<-paste0('./process/', i, '/selectKraken.reads')
  readsDat<-phylotools::read.fasta(readsPath)
  readsDat$seq.name<-gsub('\\ .*', '', readsDat$seq.name)
  lenFiltReads<-character()
  selSample<-sumDat[(sumDat$Sample.x==i),]
  targReadNames<-unique(selSample$Read)
  targReads<-readsDat[(readsDat$seq.name%in%targReadNames),]
  if (nrow(targReads) > 0) {
    for (j in 1:nrow(targReads)){
      if (nchar(targReads$seq.text[j]) > filterLength) {
        lenFiltReads<-c(lenFiltReads, targReads$seq.name[j])
      } 
    }
  }

  #result <- filter(readsDat, grepl(paste(targReadNames, collapse="|"), readsDat$seq.name))
  lenFiltReadNames<-unique(lenFiltReads)
  if (length(lenFiltReadNames) > 0) {
    outPath<-paste0('./process/', i, '/lenFiltVir.reads')
    write.table(lenFiltReadNames, file = outPath, col.names = F, row.names = F, quote = F)
  }

}

#q1<-as.data.frame(q[grepl("rs", q$SNP),])
#q2<-targReads[duplicated(targReads$seq.name)|duplicated(targReads$seq.name, fromLast=TRUE),]