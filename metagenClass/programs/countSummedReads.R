library(dplyr)
library(plyr)
library(readr)

# summary of reads per organism

setwd('/home/ubuntu/extraVol/metagenClass/2022-02-25/fastpKrakUniq/blastNtSummary/results')

# x = pattern for file, y= exact string match for virus
sumBlast<-function(x, y){
  
  
  blastNames<-c('qseqid', 'sseqid', 'stitle', 'pident', 'length', 'mismatch',
                'gapopen' , 'qstart' , 'qend', 'sstart' , 'send', 'evalue', 'bitscore')
  # final summary
  sumRes<-data.frame(matrix(ncol = 0, nrow = 0))
  # final reads table
  allReads<-data.frame(matrix(ncol = 0, nrow = 0))
  # list files
  filesList<-list.files(pattern = x)
  
  for ( i in filesList){
    blastRes<-read_delim(i, delim = "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
    resCount<-data.frame(matrix(ncol = 0, nrow = 1))
    reads<-data.frame(matrix(ncol = 0, nrow = 1))
    colnames(blastRes)<-blastNames
    blastResFilt<-blastRes[(blastRes$evalue<0.05),] # need to put evalue filtering under if statement or will cause error if no rows
    blastResFiltStr<-as.data.frame(blastResFilt[grep(y, blastResFilt$stitle, ignore.case=T),])
    if ( nrow(blastResFiltStr) > 0 ) {
      resCount$Sample<-i
      resCount$Reads_numb<-length(unique(blastResFiltStr$qseqid))
      # get unique reads
      reads<-data.frame(unique(blastResFiltStr$qseqid))
      reads$Sample<-i
      colnames(reads)[1]<-'Reads'
      reads$Virus<-y
      resCount$Virus<-y
    } else {
      resCount$Sample<-i
      resCount$Reads_numb<-0
      # get unique reads
      reads$Reads<-'No_reads'
      reads$Sample<-i
      reads$Virus<-y
      resCount$Virus<-y
    }
    sumRes<-rbind.fill(sumRes, resCount)
    allReads<-rbind.fill(allReads,  reads)
  }
  return(list(sumRes, allReads))
}

mastedSum<-sumBlast(x='mastadenovirus', y='Human mastadenovirus D')[[1]]
mastedReads<-sumBlast(x='mastadenovirus', y='Human mastadenovirus D')[[2]]

alphaSum<-sumBlast(x='alphaherpesvirus', y='human alphaherpesvirus 2')[[1]]
alphaReads<-sumBlast(x='alphaherpesvirus', y='human alphaherpesvirus 2')[[2]]

polySum<-sumBlast(x='polyomavirus', y='JC polyomavirus')[[1]]
polyReads<-sumBlast(x='polyomavirus', y='JC polyomavirus')[[2]]

#mastedSum$Virus<-'Human mastadenovirus D'
#alphaSum$Virus<-'Human alphaherpesvirus 2'
#polySum$Virus<-'JC polyomavirus'

combinedSum<-rbind(alphaSum, mastedSum, polySum)

combinedSum$Sample<-sub("_[^_]+$", "", combinedSum$Sample)

write.csv(combinedSum, './combined/blastNtSelVirResults.csv', row.names = F)


# evaluate reads 
combReads<-rbind(alphaReads, mastedReads, polyReads )
combReads$Sample<-sub("_[^_]+$", "", combReads$Sample)

write.table(combReads, './combined/blastNtSelVirReads.tsv', row.names = F, sep = '\t', quote = F)