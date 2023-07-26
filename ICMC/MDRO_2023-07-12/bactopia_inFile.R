getType1 = function(inPath, runType) {
  filesList = list.files(inPath, pattern = "_1.fastq.gz")
  if (length(filesList) > 0) {
    sample = gsub('_1.fastq.gz', "", filesList )
    r1 = normalizePath(paste0(inPath, "/",  filesList))
    r2 = normalizePath(paste0(inPath, "/", sample, "_2.fastq.gz"))
    runtype = runType
    genome_size = 0
    species = "UNKNOWN_SPECIES"
    df = data.frame(sample, runtype, genome_size, species, r1, r2)
    df$extra = ""
  }
  else {
    df = NULL
  }
  return(df)
}

getType_001 = function(inPath, runType) {
  filesList = list.files(inPath, pattern = "_R1_001.fastq.gz")
  if (length(filesList) > 0) {
    sample = gsub('_R1_001.fastq.gz', "", filesList )
    r1 = normalizePath(paste0(inPath, "/",  filesList))
    r2 = normalizePath(paste0(inPath, "/", sample, "_R2_001.fastq.gz"))
    runtype = runType
    genome_size = 0
    species = "UNKNOWN_SPECIES"
    df = data.frame(sample, runtype, genome_size, species, r1, r2)
    df$extra = ""

  } 
  else {
    df = NULL
  }
  return(df)
}

getType_R1 = function(inPath, runType) {
  filesList = list.files(inPath, pattern = "_R1.fastq.gz")
  if (length(filesList) > 0) {
    sample = gsub('_R1.fastq.gz', "", filesList )
    r1 = normalizePath(paste0(inPath, "/",  filesList))
    r2 = normalizePath(paste0(inPath, "/", sample, "_R2.fastq.gz"))
    runtype = runType
    genome_size = 0
    species = "UNKNOWN_SPECIES"
    df = data.frame(sample, runtype, genome_size, species, r1, r2)
    df$extra = ""
  } 
  else {
    df = NULL
  }
  return(df)
}

makeTables = function(curPath) {
  bactDirs = list.files(curPath)
  for (bactDir in bactDirs) {
    inPath = paste0(curPath, "/", bactDir, "/input/")
    type_1 = getType1(inPath = inPath, runType = "paired-end")
    type_001 = getType_001(inPath = inPath, runType = "paired-end")
    type_R1 = getType_R1(inPath = inPath, runType = "paired-end")
    combDat = rbind(type_1, type_001, type_R1)
    outPath = paste0(curPath, "/", bactDir, "/fastqs.txt")
    write.table(combDat, outPath, row.names = F, col.names = T, quote = F, sep = "\t")
    if ( !is.null(combDat) )  {
      if ( (length(list.files(inPath)) / 2) != nrow(combDat) ) {
        print(paste0(bactDir, "not all of the samples are recorded"))
      }
    }
  }
}

makeTables(curPath = "bactopia")