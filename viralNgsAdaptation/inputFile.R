setwd("./")

fastq_1<-list.files('./input', pattern = '_R1_')


sample<-data.frame(gsub("\\_.*", "", fastq_1))

write.table(sample, 'samples.txt', col.names = F, row.names = F, quote = F)