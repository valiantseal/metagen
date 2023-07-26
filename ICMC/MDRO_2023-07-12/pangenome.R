

prepFiles = function(inPath, exclList, exclude) {
  if (exclude == T) {
    exclSamples = read.table(exclList, F)
    exclSamples = exclSamples$V1
    exclSamples = c(exclSamples, 'bactopia-runs')
  } else {
    exclSamples = 'bactopia-runs'
  }

  bactDirs = list.files(inPath)
  for (bactDir in bactDirs) {
    newBactDir = paste0("pangenome/", bactDir, "/")
    bactOutDir = paste0("pangenome/", bactDir, "/bactopia_output/")
    
    if (!file.exists(newBactDir)) {
      dir.create(bactOutDir, recursive = T, showWarnings = F)
    }
    
    curSamples = list.files(paste0(inPath, "/", bactDir, "/output/"))
    curSamples = curSamples[!curSamples%in%exclSamples]
    
    for (curSample in curSamples) {
      targDir = paste0(inPath, "/", bactDir, "/output/", curSample, "/")
      cmd_str = paste0("mv ", targDir, " ", bactOutDir, curSample)
      try({
        system(cmd_str)
      })
    }
    
  }
}

# prepFiles(inPath = "bactopia", exclList = "dereplicatedSamples.txt", exclude = T)
# prepFiles(inPath = "bactopia_fasta", exclList = "dereplicatedSamples.txt", exclude = T)
# prepFiles(inPath = "bactopia_drep", exclList = "dereplicatedSamples.txt", exclude = F)

runPangenome = function(inPath) {
  bactDirs = list.files(inPath)
  for (bactDir in bactDirs) {
    runDir = paste0(inPath, "/", bactDir, "/")
    setwd(runDir)
    cmd_str = paste0("bactopia --wf pangenome --bactopia ./bactopia_output --outdir ./panOut --max_cpus 32 --skip_phylogeny --skip_recombination --run_name ", bactDir)
    try({
      system(cmd_str)
    })
    setwd("../../")
  }
}

runPangenome(inPath = "pangenome")