samples<-read.table('newdir.txt', F)

colnames(samples)<-'uuid'

samples$Species_full<-'Klebsiella pneumoniae'

write.csv(samples, "./metadata/metadata.csv", row.names = F)