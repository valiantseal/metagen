samplesList<-read.table('samples.list')

samplesList<-samplesList$V1

filesList<-list.files('./input', pattern = '_1.fastq.gz')

filesList<-gsub('\\_.*', '', filesList)

missingSamples<-samplesList[!(samplesList%in%filesList)]

print(paste0(length(missingSamples), ' samples missing'))

print(missingSamples)

      