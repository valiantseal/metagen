bacteria.list<-read.table('bacteria.list')

bacteria.list<-bacteria.list$V1

mlstSummary<-data.frame(matrix(ncol = 0, nrow=0))

for ( bacteria in bacteria.list){
  mlstPath<-paste0("bactopia_gtdbtk/", bacteria, "/mlstOut/bactopia-tools/mlst/mlst/mlst.tsv")
  curMlst<-read.table(mlstPath, F, sep = '\t')
  selMlst<-curMlst[, 1:3]
  mlstSummary<-rbind(mlstSummary, selMlst)
}

colnames(mlstSummary)<-c("Sample", "PubMLST", "ST")

mlstSummary$Sample<-gsub('.fna.gz', '', mlstSummary$Sample)

dir.create('custom_output', showWarnings = FALSE)

write.csv(mlstSummary, './custom_output/mlstSummary.csv', row.names = F)