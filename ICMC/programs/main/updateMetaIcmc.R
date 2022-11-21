# library(aws.s3)

df<-read.csv('./metadata/idSummary.csv')

client<-read.table('./run_info/client.name')
client<-client$V1

project<-read.table('./run_info/run.name')
project<-project$V1

outPath<-read.table('./run_info/run.path')
outPath<-outPath$V1

runDate<-file.info('download.me')$ctime

runDate<-gsub('\\ .*', '', runDate)

#datNames<-c("Sample", "Name", "Project", "Date_processed", "Gtdb_species", "Original_species", "Path")

colnames(df)[1]<-'Sample'

df$Name<-client
df$Project<-project
df$Path<-outPath
df$Date_processed<-runDate

write.csv(df, './metadata/updateMetaIcmc.csv', row.names = F)

system('aws s3 cp ./metadata/updateMetaIcmc.csv s3://transfer-files-emory/ICMC/updateMetaIcmc.csv')