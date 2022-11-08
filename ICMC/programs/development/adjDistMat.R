
fasta<-phylotools::read.fasta('./custom_output/assembled_genomes/4311_S50.fna')

sum(nchar(fasta$seq.text))

genDist<-read.table('/home/ubuntu/ICMC/Ahmed/euhm_cre/custom_output/pangenome/klebsiella-pneumoniae/core-genome.distance.tsv', T, sep = '\t')



nchar(coreGen$seq.text)[1]

# adjust distance table by the alignment size

# workflow
genDist<-read.table('/home/ubuntu/ICMC/Ahmed/euhm_cre/custom_output/pangenome/klebsiella-pneumoniae/core-genome.distance.tsv', T, sep = '\t')

alignment<-phylotools::read.fasta('./bactopia_gtdbtk/klebsiella-pneumoniae/panOut/bactopia-tools/pangenome/pangenome/core-genome.aln.gz')
alignLen<-nchar(alignment$seq.text)[1]

bacteriaList<-read.table('bacteria.list', F)
bacteriaList<-as.character(bacteriaList$V1)

genDist<-read.table('/home/ubuntu/ICMC/Ahmed/euhm_cre/custom_output/pangenome/klebsiella-pneumoniae/core-genome.distance.tsv', T, sep = '\t')

colnames(genDist)[1]<-'Samples'

geneDistEdit<-genDist[,2:ncol(genDist)]/alignLen

names(geneDistEdit) <- sub("^X", "", names(geneDistEdit))

colnames(geneDistEdit)[1:ncol(geneDistEdit)]<-genDist$Samples

# function for all samples

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