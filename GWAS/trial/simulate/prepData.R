library(vcfR)

tree<-ape::read.tree('/home/ubuntu/Bacgwasim/results_BacGWASim/simulations/genSim/phylogeny.nwk')

tree_dist<-ape::cophenetic.phylo(tree)

tree_dist<-cbind(rownames(tree_dist), tree_dist)

#snps<-read.delim('/home/ubuntu/Bacgwasim/results_BacGWASim/simulations/genSim/sims.vcf', skip=12, comment.char = '@', sep = '\t')

colnames(snps)[1]<-'#CHROM'

names(snps)[-1] <- sub(".*\\X", "", names(snps)[-1]) 

phen<-read.delim('/home/ubuntu/Bacgwasim/results_BacGWASim/simulations/phenSim/0/phenSim.phen', F, sep = ' ')

phen<-phen[, c(1,3)]

colnames(phen)<-c('samples', 'phenotype')

# save
write.table(phen, './input/pheno.tsv', row.names = F, col.names = T, sep = '\t', quote = F)

#write.table(snps, './input/snps.vcf', row.names = F, col.names = T, sep = '\t', quote = F)

#write.vcf(snps, './input/snps.vcf')

write.table(tree_dist, './input/distance.tsv', row.names = F, col.names = T, sep = '\t', quote = F)