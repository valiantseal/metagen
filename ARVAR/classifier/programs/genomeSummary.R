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
    #getDepth(curSample = i)
  }
  
  parallel::stopCluster(cl = cl)
  
}

editName = function(curSample, path) {
  curSample = gsub(path, "", curSample)
  curSample = gsub("\\/", "", curSample)
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

makeSummary = function(path, protocol) {
  allSamples = character()
  allCoverage = numeric()
  allDepth = numeric()
  filesList = list.files(path, full.names = T)
  for (curSample in filesList) {
    try ({
      sampleName = editName(curSample = curSample, path = path)
      inFile = paste0(curSample, "/", "coverage.tsv")
      curDf = read.delim(inFile)
      curCover = curDf$coverage
      curDepth = curDf$meandepth
      # combine
      allSamples = c(allSamples, sampleName)
      allCoverage = c(allCoverage, curCover)
      allDepth = c(allDepth, curDepth)
      
    })
  }
  combDf = data.frame(Sample = allSamples, Coverage = allCoverage, Mean_depth = allDepth)
  curNames = c(paste("Coverage", protocol, sep = "_"), paste("Mean_depth", protocol, sep = "_"))
  colnames(combDf)[2:3] =  curNames
  return(combDf)
}

metaSum = makeSummary(path = "Vivacilty_v1.0.1/process_par", protocol = "MetaSeq")
ampSum = makeSummary(path = "Vivacilty_v1.0.1/amp_process_par", protocol = "AmpSeq")

combSum = plyr::join(metaSum, ampSum, by = "Sample", type = 'left', match = "all")
combSum$Mean_depth_MetaSeq = format(as.numeric(combSum$Mean_depth_MetaSeq), scientific = F)

write.csv(combSum, "Ludy_overlapSamples_depthCover.csv", row.names = F)
system("aws s3 cp Ludy_overlapSamples_depthCover.csv s3://abombin/Vivacity/classifier/")

