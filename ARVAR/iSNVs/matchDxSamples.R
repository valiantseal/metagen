

combAll = function(filesList, inDir) {
  combDat = data.frame(matrix(nrow = 0, ncol = 0))
  for (curFile in filesList) {
    inFile = paste0(inDir, "/", curFile)
    curDf = read.table(inFile)
    curDf$Sample_id = curFile
    combDat = rbind(combDat, curDf)
  }
  return(combDat)
}

formatNames = function(samplesList) {
  newNames = character()
  samplesList = basename(samplesList)
  sampleNames = gsub("-", "_", samplesList)
  for ( curSample in sampleNames ) {
    sampNamelist = strsplit(curSample, "_")
    sampleName =  sampNamelist[[1]][1:4]
    sampleName = paste(sampleName, collapse = "_")
    sampleName = gsub('\\_S.*', "", sampleName)
    if (grepl("WATER", sampleName, ignore.case = T) == F) {
      newNames = c(newNames, sampleName)
    }
  }
  return(newNames)
}

# add dx fastqs data
addDx = function(mainDf, dxDf) {
  mainDf$Meta_dx = NA
  mainDf$Meta_dx_fastq = NA

  for ( i in 1:nrow(mainDf) ) {
    curSample = mainDf$Sample_id[i]
    if (!is.na(mainDf$Missing_Metaseq[i])) {
      curDxDf = dxDf[dxDf$Sample_id == curSample,]
      if (nrow(curDxDf) > 0) {
        curFiles = curDxDf$V6
        curFilesString = paste(curFiles, collapse = ";")
        curNames = unique(formatNames(samplesList = curFiles))
        curNamesStr = paste(curNames, collapse = ";")
        
      } else {
        curStrings = NA
        curNamesStr = NA
      }
      mainDf$Meta_dx[i] = curNamesStr
      mainDf$Meta_dx_fastq[i] = curFilesString
    }
  }
  return(mainDf)
}

dxMissing = function(df) {
  for ( i in 1:nrow(df) ) {
    metaseqStr = df$Missing_Metaseq[i]
    if (!is.na(metaseqStr)) {
      curMetaseq = strsplit(metaseqStr, ";")[[1]]
      curDx = strsplit(df$Meta_dx[i], ";")[[1]]
      curMissingMetaseq = character()
      for (curSeq in curMetaseq) {
        curSeqCh = strsplit(curSeq, "_")[[1]]
        curMainStr = curSeqCh[1:3]
        curMainStr = paste(curMainStr, collapse = ".")
        if (length(curSeqCh) > 3) {
          curLib = curSeqCh[4]
          curLib = gsub("[Ll]", "*",  curLib )
          curPattern = paste(curMainStr, curLib, sep = ".")
        } else {
          curPattern = curMainStr
        }

        if (!any(grepl(curPattern, curDx))) {
          curMissingMetaseq =  c(curMissingMetaseq, curSeq)
        }
        # stopped here recording if any of the samples are still missing after dx
        # if length curMissingMetaseq == 0 record as NA otherwise collapse vector into a string
      }
    }
  }
}

# current samples table
curDat = read.csv("SeqSamples_Eval_Depth10.csv")
metaSamples = curDat[!is.na(curDat$MetaseqNames),]
metaSamples = metaSamples[, 1:4]
rownames(metaSamples) = 1:nrow(metaSamples)

# get metaseq data from dx
filesList = list.files("dx_metaseq_paths")
combDat = combAll(filesList = filesList , inDir = "dx_metaseq_paths")
fastqs = combDat[grepl("fastq.gz", combDat$V6),]

# update main table with data from dx
metaDxUpdate = addDx(mainDf = metaSamples, dxDf = fastqs)
rownames(metaDxUpdate) = 1:nrow(metaDxUpdate)
