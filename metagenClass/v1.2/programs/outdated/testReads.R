selSamples<-mergeReads[(mergeReads$Sample.x=='Water-042022_S14_L001'),]

readsId<-data.frame(selSamples$Read)

write.table(readsId, 'testReads/waterReads.list', row.names = F, col.names = F, quote = F)


q1<-Water.042022_S14_L001__Human_alphaherpesvirus_1

q2<-q1[(q1$V1 %in%readsId$selSamples.Read),]

colnames(q2)<-c('qseqid', 'sseqid', 'stitle', 'pident', 'length', 'mismatch',
                'gapopen' , 'qstart' , 'qend', 'sstart' , 'send', 'evalue', 'bitscore')

write.csv(q2, 'testReads/waterReads.csv', row.names = F)