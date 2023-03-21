library(vcfR)
vcf <- read.vcfR( 'results_BacGWASim/simulations/genSim/sims.vcf', verbose = FALSE )

colnames(vcf@gt)[2001]

curNames<-colnames(vcf@gt)[2:2001]

newNames<-paste0('sample_', curNames)

newNames<-gsub('zero', '0', newNames)

newNames[1]

colnames(vcf@gt)[2:2001]<-newNames

# writes gz without option to write uncompressed
write.vcf(vcf, 'rename1.vcf.gz')

system('gzip -d rename1.vcf.gz')

system('mv rename1.vcf ./input_gwas/')