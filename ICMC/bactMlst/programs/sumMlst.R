dir.create('custom_output')

mlstSummary<-data.frame(matrix(ncol = 0, nrow=0))


mlstPath<-paste0("./mlstOut/bactopia-tools/mlst/mlst/mlst.tsv")
curMlst<-read.table(mlstPath, F, sep = '\t')
selMlst<-curMlst[, 1:3]
mlstSummary<-selMlst


colnames(mlstSummary)<-c("Sample", "PubMLST", "ST")

mlstSummary$Sample<-gsub('.fna.gz', '', mlstSummary$Sample)

dir.create('custom_output', showWarnings = FALSE)

curPath <- getwd()
curDir <- sub(".*\\/", "", curPath)


write.csv(mlstSummary, paste0('./custom_output/', curDir, '_mlstSummary.csv'), row.names = F)
write.csv(curMlst, paste0('./custom_output/', curDir, '_all_MLST_Info.csv'), row.names = F)