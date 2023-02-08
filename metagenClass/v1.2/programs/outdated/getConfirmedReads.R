df<-read.table('confirmedReads.list', T, sep = '\t')


viruses<-unique(df$Virus)
samples<-unique(df$Sample.x)

  
for (sample in samples){
  
  dfSample<-df[(df$Sample.x==sample),]
  samplePath<-paste0('../process/', sample, '/selectKraken.reads')
  fasta<-phylotools::read.fasta(samplePath)
  fasta$seq.name<-gsub("\\ .*", "", fasta$seq.name)
  
  
  for ( virus in viruses){
    dfVir<-dfSample[(dfSample$Virus==virus),]
    virusName<-gsub(" ", "_", virus)
    fileName<-paste0(sample, "_", virusName, '.fasta')
    virusReads<-fasta[(fasta$seq.name%in%dfVir$Read),]
    if (nrow(virusReads) > 0) {
      print(sample)
      print(virus)
      print(nrow(dfVir))
      print(nrow(virusReads))
      print(nchar(virusReads$seq.text))
      
      phylotools::dat2fasta(virusReads, fileName)
    }

  }
}

  