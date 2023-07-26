dir.create('bactopia_drep', showWarnings = F)

prepNewBact = function(inPath) {
  processedSamples = character()
  bactDirs = list.files(inPath, full.names = F)

  for (bactDir in bactDirs) {
    curPatients = list.files(paste0("drep/", bactDir), full.names =  T)
    for (patient in curPatients) {
      fastaList = list.files(patient, pattern = ".fna")
      if ( length(fastaList) > 1 ) {
        curSamples = gsub('.fna', '', fastaList)
        processedSamples = c(processedSamples, curSamples)
        newBactDir = paste0("bactopia_drep/", bactDir, "/input/")
        dir.create(newBactDir, recursive = T, showWarnings = F)
        drep_fasta_cmd = paste0("cp ", patient, "/drep/dereplicated_genomes/* ", newBactDir)
        try({
          system(drep_fasta_cmd)
        })
      }
    }
  }
  return(processedSamples)
}

#excludeSamples = prepNewBact("drep")

#write.table(excludeSamples, "dereplicatedSamples.txt", col.names = F, row.names = F)

addBactId = function(inPath) {
  bactDirs = list.files(inPath)
  for (bactDir in bactDirs) {
    bactid = gsub("_", " ", bactDir)
    bactid = stringr::str_to_sentence(bactid)
    outFile = paste0(inPath, "/", bactDir, "/bacteria.id")
    write.table(bactid, outFile, row.names = F, col.names = F)
  }
}

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
      prep_str = paste0("bactopia prepare -a .fna -p ./input > fastqs.txt")
      system(prep_str)
    }
    cmd_str = paste0("bactopia --samples fastqs.txt --species '", curBacteria, "' --outdir ./output --max_cpus 8")
    system(cmd_str)
    setwd("../../")
  }
}

addBactId(inPath = 'bactopia_drep')

runBactopia(inPath = "bactopia_drep", excludeDir = "", assembly = T)