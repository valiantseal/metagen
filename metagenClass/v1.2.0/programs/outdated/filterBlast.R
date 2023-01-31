library(readr)

blastRes<-read_delim('NtV4_blast.results', delim = "\t", escape_double = FALSE, col_names = FALSE, comment = "#", trim_ws = TRUE)

vir.list<-read.delim('../../../../virus.list', F, sep='\t')

blastNames<-c('qseqid', 'sseqid', 'stitle', 'pident', 'length', 'mismatch',
              'gapopen' , 'qstart' , 'qend', 'sstart' , 'send', 'evalue', 'bitscore')

colnames(blastRes)<-blastNames

reads<-unique(blastRes$qseqid)

filterReads<-function(x){
  combSamples<-data.frame(matrix(ncol = 0, nrow = 0))
  if (nrow(x) > 0)
    colnames(x)<-blastNames
    x<-x[!is.na(x$bitscore),]
    reads<-unique(x$qseqid)
    selVir<-as.data.frame(x[grep(paste(vir.list$V1,collapse="|"), x$stitle, ignore.case=T),])
    if (nrow(selVir) > 0) {
      reads<-unique(selVir$qseqid)
        for (read in reads){
          selDf<-selVir[(selVir$qseqid==read),]
          blastResFilt<-selDf[(selDf$evalue<0.05),]
          blastResFiltSyn<-as.data.frame(blastResFilt[!(grepl('Synthetic', blastResFilt$stitle, ignore.case=T)),])
          blastResFiltCl<-as.data.frame(blastResFiltSyn[!(grepl('clone', blastResFiltSyn$stitle, ignore.case=T)),])
          #readOrder<-blastResFiltCl[order(blastResFiltCl$bitscore, decreasing = T),]
          uniqueMatch<-unique(blastResFiltCl[, c('qseqid', 'stitle'),])
          combSamples<-rbind(combSamples, uniqueMatch)
        }
    }
    return(combSamples)
}


getUniqVirRead<-function(){
  seqVir<-data.frame(matrix(ncol = 0, nrow = 0))
  virus.list<-vir.list$V1
  for (virus in virus.list) {
    selVir<-as.data.frame(blastFiltRes[grep(virus, blastFiltRes$stitle, ignore.case=T),])
    if (nrow(selVir) > 0){
      selVir$Virus<-virus
      readVir<-unique(selVir[, c('qseqid', 'Virus')])
      seqVir<-rbind(seqVir, readVir)
    }
  } 
  return(seqVir)
}

blastFiltRes<-filterReads(blastRes)
seqVir<-getUniqVirRead()

if (nrow(seqVir) > 0){
  write.table(seqVir, 'blastResFilt.tsv', row.names = F, col.names = T, quote = T, sep = '\t')
}

print(paste0(getwd(), '______________done'))


