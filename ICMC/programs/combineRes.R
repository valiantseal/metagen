library(readr)
#setwd("/home/ubuntu/ICMC/cre_witt-cecile/bactopia/Enterobacter_cloacae/antimic_res")
setwd("./")

#genRef<-read_delim("/home/ubuntu/ICMC/cre_witt-cecile/bactopia/Enterobacter_cloacae/antimic_res/3653_S101-gene-report.txt", 
                   #delim = "\t", escape_double = FALSE, 
                  # trim_ws = TRUE)
#geneRefColName<-colnames(genRef)



## summary genes

genesList<-list.files(pattern = "-gene-report.txt")

combined<-data.frame(matrix(nrow=0, ncol=23))
#colnames(combined)[1]<-'sample_Id'
#colnames(combined)[2:23]<-geneRefColName

# combine tables:

for ( i in genesList){
  df<-read_delim(i, delim = "\t")
  #colnames(df)<-geneRefColName
  sampleId<-gsub("-gene-report.txt", "", i)
  df<-tibble::add_column(df, Sample=sampleId, .before = 1)
  combined<-rbind(combined, df)
}

# add species
species<-read.delim('../bacteria.id', F, sep = ',')
bactId<-species[1,1]

combined<-tibble::add_column(combined, Species=bactId, .before = 2)

bactIdFile<-gsub(" ", "_", bactId)
fileName<-paste0('./combined/', bactIdFile, '-res_genes_combined.tsv')

write.table(combined, file=fileName, col.names = T, row.names = F, quote = F, sep = '\t')


rm(list=ls())

## protein
genesList<-list.files(pattern = "-protein-report.txt")

combined<-data.frame(matrix(nrow=0, ncol=18))
#colnames(combined)[1]<-'sample_Id'
#colnames(combined)[2:23]<-geneRefColName

# combine tables:

for ( i in genesList){
  df<-read_delim(i, delim = "\t")
  #colnames(df)<-geneRefColName
  sampleId<-gsub("-protein-report.txt", "", i)
  df<-tibble::add_column(df, Sample=sampleId, .before = 1)
  combined<-rbind(combined, df)
}

# add species
species<-read.delim('../bacteria.id', F, sep = ',')
bactId<-species[1,1]

combined<-tibble::add_column(combined, Species=bactId, .before = 2)

bactIdFile<-gsub(" ", "_", bactId)
fileName<-paste0('./combined/', bactIdFile, '-res_proteins_combined.tsv')

write.table(combined, file=fileName, col.names = T, row.names = F, quote = F, sep = '\t')