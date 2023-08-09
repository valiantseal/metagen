fastqs = read.csv("metaseq_dx_fastqs.csv")

metaseq = read.csv("metaseq_libs_fail.csv")

exclList = c("EHC_C19_1407E_L2", "EHC_C19_1445Q_L2")
metaseqFilt = metaseq[!metaseq$combSeqName%in%exclList,]

metaSpike = metaseqFilt[(metaseqFilt$Spike.in.confirmed == "Y"),]
metaSpFail = metaseqFilt[!(metaseqFilt$Spike.in.confirmed == "Y"),]

convertFileSize = function(fastqs) {
  fastqs$Size = NA
  for ( i in 1:nrow(fastqs) ) {
    if (fastqs$V5[i] == "GB") {
      fastqs$Size[i] = fastqs$V4[i] * 1073741824
    } else if (fastqs$V5[i] == "MB") {
      fastqs$Size[i] = fastqs$V4[i] * 1048576
    } else if (fastqs$V5[i] == "KB") {
      fastqs$Size[i] = fastqs$V4[i] * 1024
    } else if (fastqs$V5[i] == "bytes") {
      fastqs$Size[i] = fastqs$V4[i]
    }
  }
  return(fastqs)
}

fastqs = convertFileSize(fastqs = fastqs)

getBiggestFiles = function(df, fastqs) {
  for ( i in 1:nrow(df)) {
    
    curFilesList = strsplit(df$combFastq[i], ";")[[1]]
    curFastqs = fastqs[fastqs$V6%in%curFilesList,]
    
  }
}
