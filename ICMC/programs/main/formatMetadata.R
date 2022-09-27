library(readxl)
library(tibble)
setwd('./metadata')

inFile<-list.files('./', pattern = 'xlsx')


metadata<-read_excel(inFile)

metadata<-tibble::add_column(metadata, Species_full=NA, .before = 5)


renameSpecies<-function(df){
  for (i in 1:nrow(df)) {
    if (df$Species[i]=='K. pneumoniae') {
      df$Species_full[i]<-'Klebsiella pneumoniae'
    } else if (df$Species[i]=='P. aeruginosa') {
      df$Species_full[i]<-'Pseudomonas aeruginosa'
    } else if (df$Species[i]=='K. aerogenes') {
      df$Species_full[i]<-'Klebsiella aerogenes'
    } else if (df$Species[i]=='C. freundii') {
      df$Species_full[i]<-'Citrobacter freundii'
    } else if (df$Species[i]=='R. ornithinolytica') {
      df$Species_full[i]<-'Raoultella ornithinolytica'
    } else if (df$Species[i]=='unknown') {
      df$Species_full[i]<-'unknown'
    } else if (df$Species[i]=='P. rettgeri') {
      df$Species_full[i]<-'Providencia rettgeri'
    }
  }
  colnames(df)[1]<-'uuid'
  df$uuid<-gsub('-', "_", df$uuid)
  return(df)
}

formatMetadata<-renameSpecies(metadata)

# custom formatting
numbMeta<-formatMetadata[is.na(formatMetadata$`Fastq Name`),]
numbMeta$uuid<-gsub('028_00', '', numbMeta$uuid)
selMeta<-numbMeta[(numbMeta$uuid >= 4311 & numbMeta$uuid <= 4453),]
selMeta$new_id<-gtdbtk.bac120.summary[1:16, 1]
selMeta$uuid<-selMeta$new_id
numMetFin<-selMeta[, 1:10]


charMeta<-formatMetadata[!is.na(formatMetadata$`Fastq Name`),]
charMeta$uuid<-charMeta$`Fastq Name`
charMeta$uuid<-gsub('-', '_', charMeta$uuid)
charMeta$uuid<-gsub(' ', '_', charMeta$uuid)
charMeta$uuid[(charMeta$uuid=='MU_8237')]<-'Mu_8237'
charMeta$uuid[(charMeta$uuid=='MU_8174')]<-'Mu_8174'

finalMeta<-rbind(numMetFin, charMeta)

write.csv(finalMeta, 'metadata.csv', row.names = F)



### make own metadata 
metadata<-data.frame(list.files('./process_par'))
metadata$Species_full<-'Escherichia coli'
colnames(metadata)[1]<-'uuid'
write.csv(metadata, './metadata/metadata.csv', row.names = F)