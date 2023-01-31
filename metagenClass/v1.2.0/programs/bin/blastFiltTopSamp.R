library(foreach)
library(doParallel)
library(readr)

useCores<-95

cores=detectCores()

vir.list<-read.delim('virus.list', F, sep='\t')

blastNames<-c('qseqid', 'sseqid', 'stitle', 'pident', 'length', 'mismatch',
              'gapopen' , 'qstart' , 'qend', 'sstart' , 'send', 'evalue', 'bitscore')


topNumbers<-c(1)


filterReads<-function(df, topNumb){
  x<-readr::read_delim(df, delim = "\t", escape_double = FALSE, col_names = FALSE, comment = "#", trim_ws = TRUE)
  combSamp<-data.frame(matrix(nrow = 0, ncol = 0))
  if (nrow(x) > 0) {
    colnames(x)<-blastNames
    x<-x[!is.na(x$bitscore),]
    blastResFiltSyn<-as.data.frame(x[!(grepl('Synthetic', x$stitle, ignore.case=T)),])
    blastResFiltCl<-as.data.frame(blastResFiltSyn[!(grepl('clone', blastResFiltSyn$stitle, ignore.case=T)),])
    reads<-unique(blastResFiltCl$qseqid)
    for ( read in reads){
      selRead<-blastResFiltCl[(blastResFiltCl$qseqid == read),]
      ordRead<-selRead[order(selRead$bitscore, decreasing = T),]
      topRes<- head(ordRead, topNumb) # stopped here
      topRes$Virus<-NA
      for (virSamp in vir.list$V1){
        topRes$Virus[grepl(virSamp, topRes$stitle, ignore.case = T)]<-virSamp
      } 
      filtRes<-topRes[!is.na(topRes$Virus),]
      combSamp<-rbind(combSamp, filtRes)
    }
  }
  if (nrow(combSamp) > 0) {
    combSampSel<-combSamp[, c('qseqid', 'Virus')]
    combSampUn<-unique(combSampSel)
    return(combSampUn)
  } else {
    print(paste0(df, '____no target viruses'))
  }
}


allDir<-list.files('./process')

for (sampleDir in allDir) {
  sampleDirPath<-paste0('./process/', sampleDir, '/')
  
  setwd(sampleDirPath)
  cl <- makeCluster(useCores, type = "FORK")
  registerDoParallel(cl)
  
  
  targDir<-'./splitSeq10K/'
  filesList<-list.files(targDir)
  targSamples<-paste0(targDir, filesList, '/NtV4_blast.results')
  sampleName<-read.table('sample.name', F)
  
  for (topNumber in topNumbers){
    blastComb<-data.frame(foreach(i=targSamples, .combine=rbind, .packages='readr') %dopar%{
      filterReads(i, topNumb = topNumber)
      print(i)
    })
    
    if (nrow(blastComb) > 0){
      blastComb$Sample<-sampleName$V1
      colnames(blastComb)[1]<-'Read'
      outFile<-paste0('blastResTop_', topNumber, '.tsv')
      
      write.table(blastComb, file = outFile, col.names = T, row.names = F, quote = T, sep = '\t') 
    }
  }
  
  parallel::stopCluster(cl = cl)
  
  print(paste0(sampleDirPath, '________done!'))
  setwd('../../')
}
