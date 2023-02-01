library(readxl)

metadata<-read_excel('/home/ubuntu/extraVol/GWAS/data/galardini2020/galardini2020_metadata.xlsx')

meta_ecol<-metadata[(metadata$Species == 'E. coli'),]

meta_ecol$`Strain ID`<-gsub('ECOR', 'ECOR-', meta_ecol$`Strain ID`)

genomes_list<-list.files('/home/ubuntu/extraVol/GWAS/data/galardini2020/genomes', pattern = '.fasta')

samples_list<-gsub('.fasta', '', genomes_list)

present_samples<-data.frame(meta_ecol[(meta_ecol$`Strain ID`%in%samples_list),]) # some ecor samples were skipped due to inconsistent labeling

present_genomes<-paste0(present_samples$Strain.ID, '.fasta')

present_genomes<-present_genomes[1:100]


for ( i in present_genomes){
  inPath<-paste0('../data/galardini2020/genomes/', i)
  newName<-gsub('.fasta', '.fna', i)
  outPath<-paste0('./fasta/', newName)
  syscom<-paste0('cp ', inPath, ' ', outPath)
  system(syscom)
}

# prepare input for bactopia
system('bactopia prepare ./fasta > fastqs.txt')


finalSamp<-gsub('.fasta', '', present_genomes)
runMeta<-meta_ecol[(meta_ecol$`Strain ID`%in%finalSamp),]

pheno<-runMeta[, c('Strain ID', 'Number of killed mice over 105')]
colnames(pheno)<-c('samples', 'mice_killed')

pheno$mice_killed<-gsub('\\+', '', pheno$mice_killed)

pheno$mice_killed<-as.numeric(as.character(pheno$mice_killed))

unique(is.na(pheno$mice_killed))


write.table(pheno, './gwas_input/pheno.tsv', row.names = F, col.names = T, quote = F, sep = '\t')

