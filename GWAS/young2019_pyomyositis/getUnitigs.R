library(foreach)
library(doParallel)

# needs pyseer environment to be active

selSamples = read.table("pangenome_samples.txt")

getFastaFiles = function(curSample, curPath) {
  dir.create("fasta_files", showWarnings = F)
  inPath = paste0(curPath, "/", curSample, "/main/assembler/", curSample, ".fna.gz")
  outPath = paste0("fasta_files/", curSample, ".fna.gz")
  file.copy(inPath, outPath)
  cmd_str = paste0("gzip -d ", outPath)
  system(cmd_str)
}

runPar = function(samplesList, curPath) {
  cores = 30
  cl = makeCluster(cores, type = "FORK")
  registerDoParallel(cl)
  
  foreach(i = samplesList) %dopar% {
    combDat = getFastaFiles(curSample = i, curPath = curPath)
  }
  stopCluster(cl)
}

#runPar(samplesList = selSamples$V1, curPath = "bactopia_output")

fastaPaths = normalizePath(list.files("fasta_files", full.names = T))

refPath = normalizePath("s_aureus.fasta")

write.table(fastaPaths, "untig_input.txt", col.names = F, row.names = F, quote = F)
write.table(refPath, "untig_ref.txt", col.names = F, row.names = F, quote = F)

runUnitig = function(inFile, refFile) {
  dir.create("pyseer_input")
  cmd_str = paste0("unitig-caller --call --refs ", refFile, " --reads ", inFile, " --pyseer --threads 28")
  system(cmd_str)
  file.rename("unitig_caller.pyseer", "pyseer_input/unitig_caller.pyseer")
}

runUnitig(inFile = "untig_input.txt", refFile = "untig_ref.txt")