library(foreach)
library(doParallel)

getCoverage = function(curSample) {
  curFile = paste0(curSample, "/", "output.bam")
  if (file.exists(curFile)) {
    outFile = paste0(curSample, "/", "depth.tsv")
    cmd_str = paste0("samtools depth -aa ", curFile, " -o ", outFile)
    system(cmd_str)
  }
}

getDepth = function(curSample) {
  curFile = paste0(curSample, "/", "output.bam")
  if (file.exists(curFile)) {
    outFile = paste0(curSample, "/", "coverage.tsv")
    cmd_str = paste0("Vivacilty_v1.0.1/programs/bin/samtools coverage --depth 0 ", curFile, " > ", outFile)
    system(cmd_str)
  }
}

runPar = function(filesList) {
  useCores = 6
  cl <- makeCluster(useCores, type = "FORK")
  registerDoParallel(cl)
  
  results = foreach(i=filesList) %dopar%{
    getCoverage(curSample = i)
    getDepth(curSample = i)
  }
  
  parallel::stopCluster(cl = cl)
  
}

editName = function(curSample) {
  exactSample = gsub("_", "-", curSample)
  sampNamelist = strsplit(exactSample, "-")
  sampleName =  sampNamelist[[1]][1:3]
  sampleName = paste(sampleName, collapse = "-")
  return(sampleName)
}

# metaseq
metaList = list.files("Vivacilty_v1.0.1/process_par", full.names = T)
runPar(filesList = metaList)

# ampseq
ampList = list.files("Vivacilty_v1.0.1/amp_process_par", full.names = T)
runPar(filesList = ampList)

makeSummary = function(filesList) {
  allSamples = character()
  allCoverage = numeric()
  allDepth = numeric()
  for (curSample in filesList) {
    try ({
      sampleName = editName(curSample = curSample)
      
      
    })
  }

}