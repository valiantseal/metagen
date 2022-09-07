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
    }
  }
  colnames(df)[1]<-'uuid'
  df$uuid<-gsub('-', "_", df$uuid)
  return(df)
}

formatMetadata<-renameSpecies(metadata)

write.csv(formatMetadata, 'metadata.csv', row.names = F)