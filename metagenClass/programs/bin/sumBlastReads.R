library(dplyr)
library(plyr)
library(readr)

# summary of reads per organism

setwd('./blastNtSummary/target_results')

# x = pattern for file, y= exact string match for virus
sumBlast<-function(x, y){
  
  # names of the blast columns that we have 
  blastNames<-c('qseqid', 'sseqid', 'stitle', 'pident', 'length', 'mismatch',
                'gapopen' , 'qstart' , 'qend', 'sstart' , 'send', 'evalue', 'bitscore')
  # final summary
  sumRes<-data.frame(matrix(ncol = 0, nrow = 0))
  # final reads table
  allReads<-data.frame(matrix(ncol = 0, nrow = 0))
  # list files with samples for each virus
  filesList<-list.files(pattern = x)
  
  for ( i in filesList){
    blastRes<-read_delim(i, delim = "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
    resCount<-data.frame(matrix(ncol = 0, nrow = 1))
    reads<-data.frame(matrix(ncol = 0, nrow = 1))
    if (nrow(blastRes) > 0){
      colnames(blastRes)<-blastNames
      blastResFilt<-blastRes[(blastRes$evalue<0.05),] 
    } else {
      blastResFilt<-blastRes
    }
    
    blastResFiltStr<-as.data.frame(blastResFilt[grep(y, blastResFilt$stitle, ignore.case=T),])
    blastResFiltSyn<-as.data.frame(blastResFiltStr[!(grepl('Synthetic', blastResFiltStr$stitle, ignore.case=T)),])
    blastResFiltCl<-as.data.frame(blastResFiltSyn[!(grepl('clone', blastResFiltSyn$stitle, ignore.case=T)),])
    if ( nrow(blastResFiltCl) > 0 ) {
      # count unique reads per virus per sample
      resCount$Sample<-i
      resCount$Reads_numb<-length(unique(blastResFiltCl$qseqid))
      # get unique reads per virus per sample
      reads<-data.frame(unique(blastResFiltCl$qseqid))
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
    # bind summary statistics per virus for all samples
    sumRes<-rbind(sumRes, resCount)
    # bind unique reads per virus for all samples
    allReads<-rbind(allReads,  reads)
  }
  return(list(sumRes, allReads))
}

# read list with target viruses
virusList<-read.table('../../virus.list', F, sep='\t')
virus.list<-virusList$V1

# get summaries and reads for each virus in the list across all samples
sumStat<-list()
for (virus in virus.list){
  pattern<-gsub(" ", "_", virus)
  sumVirus<-sumBlast(x=pattern, y=virus)
  sumStat<-c(sumStat, sumVirus)
}

# remove spaces and repeat viral names twice
virNames<-gsub(" ", "_", virus.list)
sumNames<-rep(virNames, each=2)

# make seuqence with increments by 2 
readsIndex<- seq(2, length(sumStat), by=2)

# name the results list with viral names
names(sumStat)<-sumNames

# mark reads tables
for (i in readsIndex){
  names(sumStat)[i]
  oldName<-names(sumStat)[i]
  newName<-paste0(oldName, "_reads")
  names(sumStat)[i]<-newName
}

readNames<-paste0(virNames, "_reads")

# bind all summary tables
sumAllViruses<-data.frame(matrix(ncol = 0, nrow = 0))

for (i in virNames){
  virSummary<-sumStat[[i]]
  sumAllViruses<-rbind(sumAllViruses, virSummary)
}

sumAllViruses$Sample<-gsub("\\__.*","", sumAllViruses$Sample)

write.csv(sumAllViruses, '../blastNtSelVirSummary.csv', row.names = F)

# bind all reads tables
sumAllReads<-data.frame(matrix(ncol = 0, nrow = 0))

for (i in readNames){
  readsSummary<-sumStat[[i]]
  sumAllReads<-rbind(sumAllReads, readsSummary)
}

sumAllReads$Sample<-gsub("\\__.*","", sumAllReads$Sample)


write.table(sumAllReads, '../blastNtSelVirReads.tsv', row.names = F, sep = '\t', quote = F)