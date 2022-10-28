# write functions sumAll() count the match for all viruses and sumSelVirus() count matches for specifc virus

setwd("./")

# read kraken results
kraken<-read.delim('./kraqSummary/krakenSelVirReads.tsv', T, sep = '\t')

# check that all reads identyfy only with one virus
krakCheck<-data.frame(table(kraken$Sample, kraken$Read))
krakExtrCh<-kraken[(kraken$Read=='NB551577:428:HLLYGAFX3:1:11110:20787:16405'),]

# create unqie column for merging 
kraken$read_sample_virus<-paste(kraken$Read, kraken$Sample, kraken$Virus, sep='_')
kraken$read_sample<-paste(kraken$Read, kraken$Sample, sep='_')

# read blast results
blast<-read.delim('./blastNtSummary/blastNtSelVirReads.tsv', T, sep = '\t')

# check that all blast reads identifies as only 1 virus
blastIDCheck<-data.frame(table(blast$Sample, blast$Reads))
blastExtrCheck<-blast[(blast$Reads=='NB551577:435:HMCY5AFX3:1:11203:26193:4115'),]

# create unqie column for merging
blast$read_sample_virus<-paste(blast$Read, blast$Sample, blast$Virus, sep='_')
blast$read_sample<-paste(blast$Read, blast$Sample, sep='_')

# rename some blast columns
colnames(blast)[3]<-'Blast_virus'
colnames(blast)[1]<-'Blast_read'

# see what blast reads batch kraken reads
combReads<-plyr::join(kraken, blast, by='read_sample', type='left', match='first')

# merge kraken and blast by unique column keeping only matching combinations
mergeReads<-merge(kraken, blast, by='read_sample_virus')

confRead<-unique(mergeReads[, c('Sample.x', 'Read', "Virus")])
write.table(confRead, './testReads/confirmedReads.list', row.names = F, col.names = T, sep = '\t')
# make sure that merged inputs are identical
identical(mergeReads$Virus, mergeReads$Blast_virus)

# make a summary on how many reads matched for each sample and virus
sampleVirus<-data.frame(table(mergeReads$Sample.x, mergeReads$Virus))
colnames(sampleVirus)<-c('Sample', 'Virus', 'Confirmed_reads')

write.csv(sampleVirus, 'blastKrakenConfirmedReadsTarget.csv', row.names = F)

# make sure no read was assigned to two viruses per sample
mergeCheck<-data.frame(table(mergeReads$Sample.x, mergeReads$Read))
