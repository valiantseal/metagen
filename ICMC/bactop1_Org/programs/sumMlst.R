dir.create('custom_output')

mlstSummary<-data.frame(matrix(ncol = 0, nrow=0))


mlstPath<-paste0("./mlstOut/bactopia-tools/mlst/mlst/mlst.tsv")
curMlst<-read.table(mlstPath, F, sep = '\t')
selMlst<-curMlst[, 1:3]
mlstSummary<-selMlst


colnames(mlstSummary)<-c("Sample", "PubMLST", "ST")

mlstSummary$Sample<-gsub('.fna.gz', '', mlstSummary$Sample)

dir.create('custom_output', showWarnings = FALSE)

write.csv(mlstSummary, './custom_output/mlstSummary.csv', row.names = F)