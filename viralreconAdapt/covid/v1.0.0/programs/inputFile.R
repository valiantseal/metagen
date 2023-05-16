setwd("./")

fastq_1<-list.files('../input', pattern = '_R1_')
fastq_2<-list.files("../input", pattern= "_R2_")

sample<-gsub("\\_.*", "", fastq_1)

inputDat<-data.frame(sample, fastq_1, fastq_2)

inputDat$fastq_1<-paste0('../input/', inputDat$fastq_1)

inputDat$fastq_2<-paste0('../input/', inputDat$fastq_2)

write.table(inputDat, 'input.csv', row.names = F, quote = F, sep=',')
