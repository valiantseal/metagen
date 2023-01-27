library(readr)


reads.list<-read.table('lenFiltVir.reads')
reads.list<-as.character(reads.list$V1)

dirList<-list.files('./splitSeq10K/', pattern = '.par')


blastComb<-data.frame(matrix(ncol = 0, nrow = 0))

for (dir in dirList) {
    dfPath<-paste0('./splitSeq10K/', dir, '/NtV4_blast.results')
    blastRes<-read_delim(dfPath, delim = "\t", escape_double = FALSE, col_names = FALSE, comment = "#", trim_ws = TRUE)
    if (nrow(blastRes) > 0 ){
      blastSel<-blastRes[!is.na(blastRes$X12),]
      blastReads<-blastSel[(blastSel$X1%in%reads.list),]
      if (nrow(blastReads) > 0){
        blastComb<-rbind(blastComb, blastReads)
      }
    }

}

write.table(blastComb, './virReadsFiltLen/blastResFilt.tsv', row.names = F, col.names = T, sep = '\t')