library(readr)



sumDat<-read.table('./testReads/confirmedReads.list', T)

samples<-unique(sumDat$Sample.x)


#i<-'B2E22-037A_S2_L001'

#df<-files.list[1]

selectTop<-function(x){
  allRes<-data.frame(matrix(ncol=0, nrow = 0))
  for ( i in samples){
    targPath<-paste0('./process/', i, '/virReadsFiltLen/')
    files.list<-list.files(targPath, pattern = '.par')
    for (df in files.list){
      dfPath<-paste0(targPath, df)
      blastRes<-read_delim(dfPath, delim = "\t", escape_double = FALSE, col_names = FALSE, comment = "#", trim_ws = TRUE, col_types = cols())
      if (nrow(blastRes) > 0) {
        blastTop<-head(blastRes, x)
        blastTop$Sample<-i
        allRes<-rbind(allRes, blastTop)
        
      }
    }
  }
  return(allRes)
}

topNumbers<-c(1)

for (topResult in topNumbers){
  
  
  allRes<-selectTop(x=topResult)
  
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
  
  curPath<-getwd()
  
  curDir<-sub(".*\\/", "", curPath)
  
  dir.create('./topMatch')
  
  # filter matched
  confReads<-combReadsSel[(combReadsSel$Match==TRUE),]
  
  confReadsSel<-unique(confReads[, c('Sample', 'Read', 'Virus')])
  
  # write filtered reads names
  
  outPath<-paste0('topMatch/', curDir,  '_krakBlastMatchLenFiltTop_', topResult, '.csv')
  
  write.csv(confReadsSel, outPath, row.names = F)
  
  # make summary
  sumTab<-data.frame(table(confReadsSel$Sample, confReadsSel$Virus))
  
  sumTabSel<-sumTab[(sumTab$Freq > 0),]
  
  colnames(sumTabSel)[1:2]<-c('Sample', 'Virus')
  
  sumPath<-paste0('topMatch/', curDir,  '_krakBlastMatchLenFilt_Sum_Top_', topResult, '.csv')
  
  write.csv(sumTabSel, sumPath, row.names = F)
  
  
  s3Path<-paste0('s3://transfer-files-emory/metagenClass/Dengue/', curDir, '/custom_output/topMatch/')
  
  s3CommandSum<-paste0('aws s3 cp ', sumPath, ' ', s3Path)
  
  s3CommandTab<-paste0('aws s3 cp ', outPath, ' ', s3Path)
  
  system(s3CommandSum)
  
  system(s3CommandTab)
}
