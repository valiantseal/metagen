dir.create('iqtree')




makeTree<-function(){
  resultDir<-paste0('./iqtree/')
  dir.create(resultDir)
  setwd(resultDir)
  alignPath<-paste0('./panOut/bactopia-tools/pangenome/pangenome/core-genome.aln.gz')
  iqCommand<-paste0('iqtree -s ', alignPath, ' -m TEST -nt 30 -bb 1000')
  system(iqCommand)
}


makeTree()
