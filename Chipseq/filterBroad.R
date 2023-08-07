filterAnnot = function(inFile, outBed, qval, outAnnot) {
  df = read.delim(inFile, F)
  df$qval = 10^-df$V9
  dfSel = df[df$qval < qval ,]
  sampleEdit = dfSel[, 1:6]
  sampleEdit$V6 = "+"
  write.table(sampleEdit, outBed, col.names = F, row.names = F, quote = F, sep = "\t")
  cmd_str = paste0("annotatePeaks.pl ", outBed, " mm10  > ", outAnnot)
  system(cmd_str)
  
}

filterAnnot(inFile = "macs2_poolAll_broad/KJ-C_peaks.broadPeak", 
            outBed = "macs2_poolAll_broad/KJ-C_peaks.bed", qval = 0.05, 
            outAnnot = "macs2_poolAll_broad/KJ-C.annotpeaks")

filterAnnot(inFile = "macs2_poolAll_fdr0.01_broad/KJ-C_peaks.broadPeak", 
            outBed = "macs2_poolAll_fdr0.01_broad/KJ-C_peaks.bed", qval = 0.01, 
            outAnnot = "macs2_poolAll_fdr0.01_broad/KJ-C.annotpeaks")
