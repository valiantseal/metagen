# write functions sumAll() count the match for all viruses and sumSelVirus() count matches for specifc virus

setwd("/home/ubuntu/extraVol/metagenClass/2022-02-25/fastpKrakUniq")

kraken<-read.delim('./kraqSummary/krakenSelVirReads.tsv', T, sep = '\t')
# check that all reads identyfy only with one virus
krakCheck<-data.frame(table(kraken$Sample, kraken$Read))
krakExtrCh<-kraken[(kraken$Read=='M01541:32:000000000-K7WN6:1:1101:12606:18147'),]

kraken$read_sample_virus<-paste(kraken$Read, kraken$Sample, kraken$Virus, sep='_')
kraken$read_sample<-paste(kraken$Read, kraken$Sample, sep='_')

blast<-read.delim('./blastNtSummary/blastNtSelVirReads.tsv', T, sep = '\t')
# check that all blast reads identifies as only 1 virus
blastIDCheck<-data.frame(table(blast$Sample, blast$Reads))
blastExtrCheck<-blast[(blast$Reads=='M01541:32:000000000-K7WN6:1:1103:24649:18804'),]

blast$read_sample_virus<-paste(blast$Read, blast$Sample, blast$Virus, sep='_')
blast$read_sample<-paste(blast$Read, blast$Sample, sep='_')

colnames(blast)[3]<-'Blast_virus'
colnames(blast)[1]<-'Blast_read'

combReads<-plyr::join(kraken, blast, by='read_sample', type='left', match='first')


mergeReads<-merge(kraken, blast, by='read_sample_virus')

identical(mergeReads$Virus, mergeReads$Blast_virus)

#noMatch<-mergeReads[!(mergeReads$Virus==mergeReads$Blast_virus),]

sampleVirus<-data.frame(table(mergeReads$Sample.x, mergeReads$Virus))

colnames(sampleVirus)<-c('Sample', 'Virus', 'Confirmed_reads')

write.csv(sampleVirus, 'blastKrakenConfirmedReadsTarget.csv', row.names = F)

mergeCheck<-data.frame(table(mergeReads$Sample.x, mergeReads$Read))
