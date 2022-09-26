library(plyr)
setwd('./bactopia_summary')

## genes

filesList<-list.files(pattern = ".tsv")

combSum<-data.frame(matrix(ncol = 0, nrow = 0))

for (i in filesList){
  df <-read.delim(i,T, sep='\t')
  print(ncol(df))
  bactName<-gsub('.tsv', '', i)
  df<-tibble::add_column(df, organism=bactName, .before = 2)
  
  combSum<-rbind.fill(combSum, df)
  
}

#finalSum<-combSum[ , colSums(is.na(combSum)) == 0]

write.csv(combSum, "./combined/summary_combined.csv", row.names = F)