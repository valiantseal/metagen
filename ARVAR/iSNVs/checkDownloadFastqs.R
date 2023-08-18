curDat = read.csv("SeqSamples_Eval_Depth10.csv")
# metaseq
missingSamples = unique(curDat$Missing_Metaseq)
missingSamples = missingSamples[!is.na(missingSamples)]

spreadList = function(x) {
  combList = character()
  for ( i in x ) {
    curList = strsplit(i, ";")[[1]]
    combList = c(combList, curList)
  } 
  return(combList)
}

curMisLibs = spreadList(missingSamples)

meta1  = list.files("metaseq_pass_fastqs")
meta2 = list.files("meta_fail_stats")
newMeta = c(meta1, meta2)

missingPerLib = function(missingList, curFiles) {
  combMissing = character()
  for ( i in missingList ) {
    curSeqCh = strsplit(i, "_")[[1]]
    curMainStr = curSeqCh[1:3]
    curMainStr = paste(curMainStr, collapse = ".")
    if (length(curSeqCh) > 3) {
      curLib = curSeqCh[4]
      #curLib = gsub("[Ll]", "*",  curLib )
      curPattern = paste(curMainStr, curLib, sep = ".")
    } else {
      curPattern = paste0(curMainStr, "(?![_-]L)")
    }
    
    if (!any(grepl(curPattern, curFiles, ignore.case = T, perl = T))) {
      combMissing =  c(combMissing, i)
    }
  }
  return(combMissing)
}

grepl("EHC-C19_2236B(?![_-]L)", "EHC-C19_2236B_1_R1", ignore.case = T, perl = T)
missingLibs = missingPerLib(missingList = curMisLibs, curFiles = newMeta)
length(missingLibs)


missingPerSamp = function(missingList, curFiles) {
  combMissing = character()
  for ( i in missingList ) {
    curSeqCh = strsplit(i, "_")[[1]]
    curMainStr = curSeqCh[1:3]
    curMainStr = paste(curMainStr, collapse = ".")
    
    if (!any(grepl(curMainStr, curFiles, ignore.case = T, perl = T))) {
      combMissing =  c(combMissing, i)
    }
  }
  return(combMissing)
}

missSamples = missingPerSamp(missingList = curMisLibs, curFiles = newMeta)




