sumDat<-read.table('./testReads/confirmedReads.list', T)

samples<-unique(sumDat$Sample.x)


#i<-'B2E22-037A_S2_L001'

#df<-files.list[1]

allRes<-data.frame(matrix(ncol=0, nrow = 0))

for ( i in samples){
  targPath<-paste0('./process/', i, '/virReadsFiltLen/')
  files.list<-list.files(targPath, pattern = '.par')
  for (df in files.list){
    dfPath<-paste0(targPath, df)
    blastRes<-read_delim(dfPath, delim = "\t", escape_double = FALSE,  col_names = FALSE, trim_ws = TRUE, skip = 1)
    if (nrow(blastRes) > 0){
      blastTop<-head(blastRes, 1)
      blastTop$Sample<-i
      allRes<-rbind(allRes, blastTop)
    }
  }
}

blastNames<-c('qseqid', 'sseqid', 'stitle', 'pident', 'length', 'mismatch',
                            'gapopen' , 'qstart' , 'qend', 'sstart' , 'send', 'evalue', 'bitscore')

colnames(allRes)[1:13]<-blastNames

allRes$Organism<-stringr::word(allRes$stitle, 1,2, sep=" ")

allResSel<-allRes[, c('Sample', 'qseqid', 'stitle')]

allResSel$read_sample<-paste(allResSel$qseqid, allResSel$Sample, sep = '_')

sumDat$read_sample<-paste(sumDat$Read, sumDat$Sample.x, sep = '_')

combReads<-plyr::join(allResSel, sumDat, by='read_sample', type='left', match='all')

combReadsSel<-combReads[, c('Sample', 'Read', 'stitle', 'Virus')]

for (i in 1:nrow(combReadsSel)){
  combReadsSel$Match[i]<-grepl(combReadsSel$Virus[i], combReadsSel$stitle[i], ignore.case = T)
}

write.csv(combReadsSel, 'krakBlastMatchLenFiltTop1.csv', row.names = F)