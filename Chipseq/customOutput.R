setwd("/home/flyhunter/Kai/Chipseq/pnas/")

taretDir = "custom_output/"
dir.create(taretDir, showWarnings = F)

filesList = list.files("/home/flyhunter/Kai/Chipseq/pnas/output_merge_control/bwa/mergedLibrary/macs2/narrowPeak/", pattern = "annotatePeaks.txt")

for ( i in filesList) {
  targFile = paste0("/home/flyhunter/Kai/Chipseq/pnas/output_merge_control/bwa/mergedLibrary/macs2/narrowPeak/", i)
  df = read.delim(targFile)
  print(i)
  print("Chd7" %in% df$Gene.Name)
  colnames(df)[1] = "Peaks"
  fileName = gsub(".txt", "", i)
  write.csv(df, paste0(taretDir, i, ".csv"), row.names = F)
}

filesList = list.files('/home/flyhunter/Kai/Chipseq/pnas/output_merge_control/bwa/mergedLibrary/macs2/narrowPeak/consensus')

for ( i in filesList) {
  targFile = paste0("/home/flyhunter/Kai/Chipseq/pnas/output_merge_control/bwa/mergedLibrary/macs2/narrowPeak/consensus/", i, "/", i, ".consensus_peaks.annotatePeaks.txt")
  df = read.delim(targFile)
  print(i)
  print("Chd7" %in% df$Gene.Name)
  colnames(df)[1] = "Peaks"
  fileName = paste0(i, ".consensus_peaks.annotatePeaks.csv")
  write.csv(df, file = fileName, row.names = F)
}