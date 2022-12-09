dir.create('iqtree')

bacteria.list<-read.table('bacteria.list')

bacteria.list<-bacteria.list$V1



makeTree<-function(check_samp){
  for ( i in bacteria.list){
    resultDir<-paste0('./iqtree/', i)
    dir.create(resultDir)
    setwd(resultDir)
    alignPath<-paste0('../../bactopia_gtdbtk/', i, '/panOut/bactopia-tools/pangenome/pangenome/core-genome.aln.gz')
      if (check_samp == T) {
        samples<-read.table(paste0('./test_gtdbtk/buckets/', i))
        sampNumb<-nrow(samples)
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
        
      } else {
        iqCommand<-paste0('iqtree -s ', alignPath, ' -m TEST -nt 16 -bb 1000 -pre ', i)
        system(iqCommand)
      }
    setwd('../../')
  }
}


makeTree(check_samp = F)
