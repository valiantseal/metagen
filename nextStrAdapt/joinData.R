
setwd("/home/ubuntu/strain/processData/georgia/2021-03-01_2022-06-23")

# merge metadata
metadatMerge<-function(){
  sample.metadata <- read.delim("1655840417523.metadata.tsv")
  metadatList<-list.files(pattern = "tsv")
  metadatSum<-data.frame(matrix(ncol=28))
  colnames(metadatSum)<-colnames(sample.metadata)
  for (i in metadatList){
    metadata<-read.delim(i, T, sep='\t')
    metadatSum<-rbind(metadatSum, metadata)
  }
  metaDatFilr<-metadatSum[!(is.na(metadatSum$strain)),]
  return(metaDatFilr)
}

# merge fasta
fastaMerge<-function(){
  fastaSum<-data.frame(matrix(ncol=2))
  colnames(fastaSum)<-c('seq.name', 'seq.text')
  fastaList<-list.files(pattern = '.fasta')
  for (i in fastaList){
    fasta<-phylotools::read.fasta(i)
    fastaSum<-rbind(fastaSum, fasta)
  }
  fastaFiltr<-fastaSum[!(is.na(fastaSum$seq.name)),]
  return(fastaFiltr)
}

# remove duplicated samples
metadata<-unique(metadatMerge())
metadataUniq<-metadata[!duplicated(metadata$strain), ]
fasta<-unique(fastaMerge())


setwd('./joined')

write.table(metadataUniq, 'georgia_2021-03-01_2022-06-23.tsv', row.names = F, quote = F, sep = '\t')

phylotools::dat2fasta(fasta, 'georgia_2021-03-01_2022-06-23.fasta')