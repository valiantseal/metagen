makeTree = function(check_samp, inPath) {
  bactDirs = list.files(inPath)
  for (bactDir in bactDirs) {
    curDirs = list.files(paste0(inPath, "/", bactDir, "/panOut/bactopia-runs"), full.names = T)
    # change next line when run with normal bactopia runs
    samples = list.files(paste0(inPath, "/", bactDir, "/bactopia_output/"))
    samples = samples[!samples%in%'bactopia-runs']
    samplesNumb = length(samples)
    for (curDir in curDirs) {
      setwd(curDir)
      if (check_samp == T) {
        if (samplesNumb > 3) {
          iqCommand<-paste0('iqtree -s core-genome.aln.gz -m TEST -nt 16 -bb 1000')
        } else if (samplesNumb == 3) {
          iqCommand<-paste0('iqtree -s core-genome.aln.gz -m TEST -nt 16')
        }
        try({
          system(iqCommand)
        })
      } else {
        iqCommand<-paste0('iqtree -s core-genome.aln.gz -m TEST -nt 16 -bb 1000')
        try({
          system(iqCommand)
        })
      }
      
      setwd("../../../../../")
    }
  }
}

makeTree(check_samp = T, inPath = "pangenome")