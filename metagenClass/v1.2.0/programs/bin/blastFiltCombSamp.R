library(foreach)
library(doParallel)
library(readr)

useCores<-95

cores=detectCores()

vir.list<-read.delim('virus.list', F, sep='\t')

blastNames<-c('qseqid', 'sseqid', 'stitle', 'pident', 'length', 'mismatch',
              'gapopen' , 'qstart' , 'qend', 'sstart' , 'send', 'evalue', 'bitscore')




filterReads<-function(df){
  x<-readr::read_delim(df, delim = "\t", escape_double = FALSE, col_names = FALSE, comment = "#", trim_ws = TRUE)
  combSamples<-data.frame(matrix(ncol = 0, nrow = 0))
  if (nrow(x) > 0)
    colnames(x)<-blastNames
  x<-x[!is.na(x$bitscore),]
  reads<-unique(x$qseqid)
  selVir<-as.data.frame(x[grep(paste(vir.list$V1,collapse="|"), x$stitle, ignore.case=T),])
  if (nrow(selVir) > 0) {
      blastResFilt<-selVir[(selVir$evalue<0.05),]
      blastResFiltSyn<-as.data.frame(blastResFilt[!(grepl('Synthetic', blastResFilt$stitle, ignore.case=T)),])
      blastResFiltCl<-as.data.frame(blastResFiltSyn[!(grepl('clone', blastResFiltSyn$stitle, ignore.case=T)),])
      #readOrder<-blastResFiltCl[order(blastResFiltCl$bitscore, decreasing = T),]
      
      # code to add formal virus name
      blastResFiltCl$Virus<-NA
      for (virSamp in vir.list$V1){
        blastResFiltCl$Virus[grepl(virSamp, blastResFiltCl$stitle, ignore.case = T)]<-virSamp
      }
      uniqueMatch<-blastResFiltCl[, c('qseqid', 'Virus')]
      combSamples<-rbind(combSamples, uniqueMatch)
  }
  combSamplesUn<-unique(combSamples)
  if (nrow(combSamplesUn) > 0) {
    return(combSamplesUn)
  } else {
    print(' no target viruses')
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
  
  blastComb<-data.frame(foreach(i=targSamples, .combine=rbind, .packages='readr') %dopar%{
    filterReads(i)
  })
  
  sampleName<-read.table('sample.name', F)
  
  blastComb$Sample<-sampleName$V1
  
  write.table(blastComb, 'blastTargVirRes.tsv', col.names = T, row.names = F, quote = T, sep = '\t')
  
  
  parallel::stopCluster(cl = cl)
  
  print(paste0(sampleDirPath, '________done!'))
  setwd('../../')
}

