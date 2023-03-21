phen<-read.delim('results_BacGWASim/simulations/phenSim/0/phenSim.phen', F, sep = ' ')

phen<-phen[, c(1,3)]

colnames(phen)<-c('samples', 'phenotype')

# save
write.table(phen, 'pheno_0.tsv', row.names = F, col.names = T, sep = '\t', quote = F)