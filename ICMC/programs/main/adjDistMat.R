bacteriaList<-read.table('bacteria.list', F)
bacteriaList<-as.character(bacteriaList$V1)


adjDistAlignLen<-function(x){
  for (i in x) {
    genDistPath<-paste0('./custom_output/pangenome/', i, '/core-genome.distance.tsv')
    if (file.exists(genDistPath)) {
      genDist<-read.table(genDistPath, T, sep = '\t')
      colnames(genDist)[1]<-'Samples'
      names(genDist) <- sub("^X", "", names(genDist))
      alignmentPath<-paste0('./bactopia_gtdbtk/', i, '/panOut/bactopia-tools/pangenome/pangenome/core-genome.aln.gz')
      alignment<-phylotools::read.fasta(alignmentPath)
      alignLen<-nchar(alignment$seq.text)[1]
      geneDistEdit<-genDist[,2:ncol(genDist)]/alignLen
      colnames(geneDistEdit)[1:ncol(geneDistEdit)]<-genDist$Samples
      adjDist<-cbind(genDist$Samples, geneDistEdit)
      colnames(adjDist)<-colnames(genDist)
      
      adjDistPath<-paste0('./custom_output/pangenome/', i, '/core-genome_distance_adj_alignLen.tsv')
      write.table(adjDist, adjDistPath, row.names = F, col.names = T, sep = '\t')
      
    }
  }
}

adjDistAlignLen(bacteriaList)

