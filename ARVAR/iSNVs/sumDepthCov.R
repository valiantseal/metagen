library(doParallel)
library(foreach)


getFilesList = function(inPath, pattern) {
  filesList = list.files(inPath, pattern = pattern, full.names = T)
  return(filesList)
}

ampList = getFilesList(inPath = "/home/ubuntu/extraVol/Viralrecon/covid/iSNVs_ampseq_2023-08-22/output/variants/bowtie2", pattern = "_trim.sorted.bam$")
metaList = getFilesList(inPath = "/home/ubuntu/extraVol/Viralrecon/covid/iSNVs_metaseq_2023-08-22/output/variants/bowtie2", pattern = ".sorted.bam$")

getStats = function(inFile) {
  cmd_str = paste0("/home/ubuntu/extraVol/bio-soft/samtools coverage --depth 0 ", inFile)
  statsStr = system(cmd_str, intern = T)
  outStrs = strsplit(statsStr, " ")
  curStats = strsplit(outStrs[[2]], "\t")[[1]]
  curCov = curStats[6]
  curDepth = curStats[7]
  curFile = basename(inFile)
  curFile = gsub("\\..*", "", curFile)
  df = data.frame(curFile, curCov, curDepth)
  colnames(df) = c("Sample_id", "Coverage", "Mean_depth")
  return(df)
}

# run parallel
runPar = function(filesList) {
  useCores = 6
  cl <- makeCluster(useCores, type = "FORK")
  registerDoParallel(cl)
  
  combDat = foreach(i = filesList, .combine  = 'rbind') %dopar% {
    getStats(inFile = i)
  }
  parallel::stopCluster(cl = cl)
  return(combDat)
}

# ampStats = runPar(filesList = ampList)
# write.csv(ampStats, "ampseqCovDepth.csv", row.names = F)

metaStats = runPar(filesList = metaList)
write.csv(metaStats, "metaseqCovDepth.csv", row.names = F)