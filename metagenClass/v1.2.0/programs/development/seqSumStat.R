
#system("seqkit stats -a ./process/B2E22-042A_S5_L001/merged_reads.fa > seqkit.stats")

samples.list<-read.table('newdir.list')
samples.list<-samples.list$V1

combDf<-data.frame(matrix(ncol = 0, nrow = 0))

for (i in samples.list){
  seqPath<-paste0('./process/', i, '/merged_reads.fa')
  seqCommand<-paste0('seqkit stats -a ', seqPath, ' > seqkit.stats')
  system(seqCommand)
  seq.result<-read.table('seqkit.stats', T)
  combDf<-rbind(combDf, seq.result)
}

median(combDf$Q2)


prevMed<-median(combDf$Q2)