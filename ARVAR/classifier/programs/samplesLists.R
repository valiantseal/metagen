# overlap
overlapSamplesDf = read.csv("result_tables/Ludy_metaAmpIvar_overlapSnv.csv")
overlapSamples = unique(overlapSamplesDf$Sample)

# metaseq
metaseq1 = read.csv("/home/ubuntu/extraVol/Viralrecon/covid/Ludy_metaseq_2023-05-16/input.csv")
metaseq2 = read.csv("/home/ubuntu/extraVol/Viralrecon/covid/Ludy_metaseq_2023-05-24/input.csv")
metaComb = unique(c(metaseq1$sample, metaseq2$sample))

#ampSeq 
amp1 = read.csv("/home/ubuntu/extraVol/Viralrecon/covid/Ludy_Apr242023/input.csv")
amp2 = read.csv("/home/ubuntu/extraVol/Viralrecon/covid/Ludy_ampseq_2023-06-16/input.csv")
ampComb = unique(c(amp1$sample, amp2$sample))

write.csv(metaComb, "Ludy_metaseq_samples.csv", row.names = F)
write.csv(ampComb, "Ludy_ampseq_samples.csv", row.names = F)

system("aws s3 cp Ludy_metaseq_samples.csv s3://abombin/Vivacity/classifier/")
system("aws s3 cp Ludy_ampseq_samples.csv s3://abombin/Vivacity/classifier/")
###



editNames = function(x) {
  newNames = character()
  listEd = gsub("_", "-", x)
  listEd = toupper(listEd)
  combSamp = strsplit(listEd, "-")
  for ( i in 1:length(combSamp)) {
    curName = paste(combSamp[[i]][1:3], collapse = "-")
    newNames = c(newNames, curName)
  }
  return(newNames)
}

ampCombEd = editNames(ampComb)
metaCombEd = editNames(metaComb)

newOverlap = metaCombEd[(metaCombEd%in%ampCombEd)]

missOverlap = newOverlap[!(newOverlap%in%overlapSamples)]

nonOverlapMeta = metaCombEd[!(metaCombEd%in%ampCombEd)]

# write edited samples names to folder where I will run viralrecon

write.csv(metaCombEd, "/home/ubuntu/extraVol/Viralrecon/covid/Ludy_metaseq_2023-06-27/Ludy_metaseq_samples.csv", row.names = F)
write.csv(ampCombEd, "/home/ubuntu/extraVol/Viralrecon/covid/Ludy_metaseq_2023-06-27/Ludy_ampseq_samples.csv", row.names = F)