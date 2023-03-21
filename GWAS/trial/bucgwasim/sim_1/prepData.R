dir.create('input_gwas')


# rename alignment
curNames<-0:1999

newNames<-paste0('sample_',curNames)

renameDf<-data.frame(cbind(curNames, newNames))

system.time({phylotools::rename.fasta(infile = 'results_BacGWASim/simulations/genSim/genSim.fasta', 
                         ref_table=renameDf, outfile = "alignment.fasta")})

# make a reference sequence for unitigs

fasta<-phylotools::read.fasta('results_BacGWASim/simulations/genSim/genSim.fasta')

reference<-fasta[1,]

reference[1,1]<-'reference_seq'

phylotools::dat2fasta(reference, 'reference.fasta')

# rename tree tips
tree<-ape::read.tree('results_BacGWASim/simulations/genSim/phylogeny.nwk')

curTips<-tree$tip.label

newTips<-paste0('sample_', curTips)

tree$tip.label<-newTips

tree$tip.label<-gsub('zero', '0', tree$tip.label)

plot(tree)

ape::write.tree(tree, 'tree_rename.nwk')

# edit phenotype data

pheno<-read.table('results_BacGWASim/simulations/phenSim/0/phenSim.phen', F)

pheno<-pheno[, 2:3]

colnames(pheno)<-c('samples', 'phenotype')

identical(pheno$samples[2:1999], renameDf$curNames[2:1999])

pheno$samples<-renameDf$newNames

write.table(pheno, './input_gwas/phnotype_0.tsv', col.names = T, row.names = F, quote = F, sep = '\t')

# split alignment in separate fasta for unitig caller

fasta$seq.name<-newNames

dir.create('split_fasta')


for ( i in newNames ) {
  curSeq<-fasta[(fasta$seq.name == i ),]
  outFile<-paste0('split_fasta/' ,i, '.fasta')
  phylotools::dat2fasta(curSeq, outfile = outFile)
}

# make unitig caller input files
curPath<-getwd()
genomFiles<-list.files('split_fasta', full.names = T)
genomFiles<-paste(curPath, genomFiles, sep = '/')
refPath<-paste(curPath, 'reference.fasta', sep = '/')

write.table(genomFiles, 'untig_input.txt', col.names = F, row.names = F, quote = F)

write.table(refPath, 'untig_ref.txt', col.names = F, row.names = F, quote = F)