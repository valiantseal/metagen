library(dplyr)

fastqs = read.csv("ampseq_combined_paths.csv")
fastqs$Pbnas_path = gsub('\\\\', "/", fastqs$Pbnas_path)
fastqs = fastqs[grepl("fastq.gz", fastqs$Pbnas_path),]

ampseq = read.csv("ampseq_libs_fail.csv")


ampSpike = ampseq[(ampseq$Spike.in.confirmed == "AMPSEQ"),]
ampSpFail = ampseq[!(ampseq$Spike.in.confirmed == "AMPSEQ"),]
rownames(ampSpFail) = NULL

getBiggestFiles = function(df, fastqs) {
  combDat = data.frame(matrix(nrow = 0, ncol = 0))
  for ( i in 1:nrow(df)) {
    curFilesList = strsplit(df$combFastq[i], ";")[[1]]
    curFastqs = fastqs[fastqs$Pbnas_path%in%curFilesList,]
    curFastqs$Group = NA
    for (j in 1:nrow(curFastqs)) {
      pathList = strsplit(curFastqs$Pbnas_path[j], "/")[[1]]
      grInd = length(pathList) - 1
      elem =  pathList[grInd]
      curFastqs$Group[j] =  elem
    }

    
    size_mean <- curFastqs %>%
      group_by(Group) %>%
      summarise(SumValue = mean(File_size))
    max_mean = max(size_mean$SumValue)
    max_df = size_mean[size_mean$SumValue ==  max_mean,]
    max_gr = max_df$Group[1]
    cur_fastqs_max = curFastqs[curFastqs$Group ==  max_gr,]
    cur_fastqs_max$Lib_name = df$combSeqName[i]
    combDat = rbind(combDat, cur_fastqs_max)
  }
  return(combDat)
}

bigamp = getBiggestFiles(df = ampSpike, fastqs = fastqs)
bmSum = data.frame(table(bigamp$Lib_name))

bigampFail = unique(getBiggestFiles(df = ampSpFail, fastqs = fastqs))


splitFailFiles = function(df, fastqs) {
  combDat = data.frame(matrix(nrow = 0, ncol = 0))
  for ( i in 1:nrow(df)) {
    curFilesList = strsplit(df$combFastq[i], ";")[[1]]
    curFastqs = fastqs[fastqs$Pbnas_path%in%curFilesList,]
    curFastqs$Group = NA
    for (j in 1:nrow(curFastqs)) {
      pathList = strsplit(curFastqs$Pbnas_path[j], "/")[[1]]
      grInd = length(pathList) - 1
      elem =  pathList[grInd]
      curFastqs$Group[j] =  elem
    }
    combDat = rbind(combDat, curFastqs)
  }
  return(combDat)
}

failedFilesFastqs = unique(splitFailFiles(df = ampSpFail, fastqs = fastqs))
rownames(failedFilesFastqs) = NULL

write.csv(failedFilesFastqs, "ampseq_failed_diffSize.csv", row.names = F)
system("aws s3 cp ampseq_failed_diffSize.csv s3://abombin/ARVAR/iSNVs/")