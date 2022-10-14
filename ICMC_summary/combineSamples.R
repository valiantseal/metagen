setwd('/home/ubuntu/ICMC/summaries')

filesList<-list.files(pattern='.csv')

allSamples<-data.frame(matrix(ncol=0, nrow = 0))

for (i in filesList){
  df<-read.csv(i)
  allSamples<-rbind(df, allSamples)
}

allSamples$Date_processed<-gsub("\\ .*", "", allSamples$Date_processed)

allSamples$Gtdb_species[allSamples$Gtdb_species=='unkown' | allSamples$Gtdb_species=='unknown']<-NA

length(unique(allSamples$Sample))

nameProjSum<-data.frame(table(allSamples$Name, allSamples$Project))

nameProjFilt<-nameProjSum[!(nameProjSum$Freq==0),]

colnames(nameProjFilt)<-c('Name', 'Project', 'Samples_Number')

write.csv(allSamples, './combined/allSamples.csv', row.names = F)

write.csv(nameProjFilt, './combined/allSamplesNameProjSum.csv', row.names = F)