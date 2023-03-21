filesList<-list.files('results_BacGWASim/simulations/phenSim/')

dir.create('./phenotypes/')

for ( i in filesList) {
  targPhen <- paste0('results_BacGWASim/simulations/phenSim/', i, '/phenSim.phen')
  phen<-read.delim(targPhen, F, sep = ' ')
  phen<-phen[, c(1,3)]
  colnames(phen)<-c('samples', 'phenotype')
  outFile<-paste0('./phenotypes/', 'pheno_', i, '.tsv')
  write.table(phen, outFile, row.names = F, col.names = T, sep = '\t', quote = F)
}


