runBactopia = function(inPath, excludeDir, assembly) {
  dirsList = list.files(inPath, full.names = T)
  filtDirs = dirsList[!(dirsList%in%excludeDir)]
    for (curDir in dirsList) {
      setwd(curDir)
      try({
        system("rm -rf output work")
      })
      curBacteria = read.table("bacteria.id", sep = "\t")
      curBacteria = curBacteria$V1
      if (assembly == T) {
        prep_str = paste0("bactopia prepare -a .fasta -p ./input > fastqs.txt")
        system(prep_str)
      }
      cmd_str = paste0("bactopia --samples fastqs.txt --species '", curBacteria, "' --outdir ./output --max_cpus 8")
      system(cmd_str)
      setwd("../../")
  }
}


#runBactopia(inPath = "bactopia", excludeDir = "", assembly = F)

runBactopia(inPath = "bactopia_fasta", excludeDir = "", assembly = T)

