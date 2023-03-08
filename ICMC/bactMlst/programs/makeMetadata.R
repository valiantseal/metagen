library(stringr)

curDate = Sys.Date()

curPath = getwd()

curDir = sub(".*\\/", "", curPath)

dirSplit = strsplit( curPath , "\\/" )
curName = dirSplit[[1]]
curName = curName[length(curName)-1]

samples = list.files('./input/', pattern = '_R1.fastq.gz')
samples = gsub('_R1.fastq.gz', '', samples)

gtdb_species = 'Not_tested'
original_species = 'Haemophilus influenzae'

curPath = 's3://transfer-files-emory/ICMC/Sarah/hflu_2023-02-27'

metadata = data.frame(samples)

metadata$Name = curName
metadata$Project = curDir
metadata$Date_processed = curDate
metadata$Gtdb_species = gtdb_species
metadata$Original_species = original_species
metadata$Path = curPath
colnames(metadata)[1] = 'Sample'

dir.create('metadata')

write.csv(metadata, './metadata/metadata.csv', row.names = F)

write.csv(metadata, './metadata/updateMetaIcmc.csv', row.names = F)

system('aws s3 cp ./metadata/updateMetaIcmc.csv s3://transfer-files-emory/ICMC/updateMetaIcmc.csv')

write.table(samples, 'samples.list', row.names = F, col.names = F, quote = F)