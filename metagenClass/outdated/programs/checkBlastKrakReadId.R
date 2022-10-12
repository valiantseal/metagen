# write functions sumAll() count the match for all viruses and sumSelVirus() count matches for specifc virus


kraken<-read.delim('/home/ubuntu/extraVol/metagenClass/2022-02-25/fastpKrakUniq/kraqSummary/combined/krakenSelVirReads.tsv', T, sep = '\t')
kraken$read_sample_virus<-paste(kraken$Read, kraken$Sample, kraken$Virus, sep='_')

blast<-read.delim('/home/ubuntu/extraVol/metagenClass/2022-02-25/fastpKrakUniq/blastNtSummary/results/combined/blastNtSelVirReads.tsv', 
                  T, sep = '\t')

blast$read_sample_virus<-paste(blast$Read, blast$Sample, blast$Virus, sep='_')

colnames(blast)[3]<-'Blast_virus'
colnames(blast)[1]<-'Blast_read'

combReads<-plyr::join(kraken, blast, by='read_sample_virus', type='left', match='first')


mergeReads<-merge(kraken, blast, by='read_sample_virus')

identical(mergeReads$Virus, mergeReads$Blast_virus)

#noMatch<-mergeReads[!(mergeReads$Virus==mergeReads$Blast_virus),]

sampleVirus<-data.frame(table(mergeReads$Sample.x, mergeReads$Virus))


# make sure results match original blast results

