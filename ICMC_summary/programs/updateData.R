
allSamples<-read.csv('../combined/allSamples.csv')

system('aws s3 cp s3://transfer-files-emory/ICMC/updateMetaIcmc.csv ../combined/')


updateDf<-read.csv('../combined/updateMetaIcmc.csv')

updMeta<-rbind(allSamples, updateDf)

write.csv(updMeta, '../combined/allSamples.csv', row.names = F)