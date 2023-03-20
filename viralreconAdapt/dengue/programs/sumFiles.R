mosList = list.files('./summary', pattern ='summary', full.names = T)

mosSum = data.frame(matrix(nrow = 0, ncol = 0))

for (i in mosList){
  mos = read.delim(i, T)
  curFile = gsub('./summary/', '', i)
  curFile = gsub('.mosdepth.summary.txt', '', curFile)
  mos$File = curFile
  mosSum = rbind(mosSum, mos)
}

mosFilt = mosSum[(mosSum$chrom == 'total'),]

write.csv(mosFilt, 'hrtv_2023-03-23_depthSummary.csv', row.names = F)

# coverage 
covList = list.files('./summary', pattern ='coverage', full.names = T)

covSum = data.frame(matrix(nrow = 0, ncol = 0))

for (i in covList){
  cov = read.delim(i, T)
  curFile = gsub('./summary/', '', i)
  curFile = gsub('.coverage.tsv', '', curFile)
  cov$File = curFile
  covSum = rbind(covSum , cov)
}

write.csv(covSum, 'hrtv_2023-03-23_coverageSummary.csv', row.names = F)

