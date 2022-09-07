library(readr)

## genes

filesList<-list.files(pattern = "genes_combined.tsv")

combGene<-data.frame(matrix(ncol=24, nrow=0))

for (i in filesList){
  df <-read_delim(i, delim='\t')
  combGene<-rbind(combGene, df)
  
}

## proteins

protList<-list.files(pattern = "proteins_combined.tsv")

combProt<-data.frame(matrix(ncol=20, nrow=0))

for (j in protList){
  data <-read_delim(j, delim='\t')
  combProt<-rbind(combProt, data)
  
}

write.csv(combGene, "./combined/antimicRes_genes_combined.csv", row.names = F)


write.csv(combProt, "./combined/antimicRes_proteins_combined.csv", row.names = F)
