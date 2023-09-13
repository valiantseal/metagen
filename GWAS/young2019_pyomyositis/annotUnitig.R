library(foreach)
library(doParallel)

# needs pyseer environment to be active

selSamples = read.table("pangenome_samples.txt")

getGffFiles = function(curSample, curPath) {
  dir.create("gff_files", showWarnings = F)
  inPath = paste0(curPath, "/", curSample, "/main/annotator/prokka/", curSample, ".gff.gz")
  outPath = paste0("gff_files/", curSample, ".gff.gz")
  file.copy(inPath, outPath)
  cmd_str = paste0("gzip -d ", outPath)
  system(cmd_str)
}

runPar = function(samplesList, curPath) {
  cores = 15
  cl = makeCluster(cores, type = "FORK")
  registerDoParallel(cl)
  
  foreach(i = samplesList) %dopar% {
    combDat = getGffFiles(curSample = i, curPath = curPath)
  }
  stopCluster(cl)
}

#runPar(samplesList = selSamples$V1, curPath = "bactopia_output")

##
makeRefFile = function() {
  fastaList = character()
  gffList = character()
  for (curSample in selSamples$V1) {
    curFasta = list.files("fasta_files/", pattern = curSample)
    curGff = list.files("gff_files/", pattern = curSample)
    fastaList = c(fastaList, curFasta)
    gffList = c(gffList, curGff)
  }
  combDat = data.frame(fastaList, gffList)
  colnames(combDat) = c("V1", "V2")
  combDat$V3 = "draft"
  baseRef = read.delim("annotate/references.txt", F)
  finalDat = rbind(baseRef, combDat)
  return(finalDat)
}

refDf = makeRefFile()

