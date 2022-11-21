dir.create('iqtree')

bacteria.list<-read.table('bacteria.list')

bacteria.list<-bacteria.list$V1

for ( i in bacteria.list){
  resultDir<-paste0('./iqtree/', i)
  dir.create(resultDir)
  # check number of samples
  samples<-read.table(paste0('./test_gtdbtk/buckets/', i))
  sampNumb<-nrow(samples)
  setwd(resultDir)
  alignPath<-paste0('../../bactopia_gtdbtk/', i, '/panOut/bactopia-tools/pangenome/pangenome/core-genome.aln.gz')
  if (sampNumb == 3 ) {
    iqCommand<-paste0('iqtree -s ', alignPath, ' -m TEST -nt 16 -pre ', i)
    system(iqCommand)
  } else if (sampNumb > 3 ){
    iqCommand<-paste0('iqtree -s ', alignPath, ' -m TEST -nt 16 -bb 1000 -pre ', i)
    system(iqCommand)
  } else {
    fileConn<-file("output.txt")
    writeLines(c("Less than 3 samples"), fileConn)
    close(fileConn)
  }
  setwd('../../')
}
