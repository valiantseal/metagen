
makeBed = function(inFile, outFile, header) {
  df = read.delim(inFile, header)
  dfSub = df[, 1:3]
  colnames(dfSub) = c("chromosome", "start", "end")
  dfSub$peak_id = NA
  for ( i in 1:nrow(dfSub) ) {
    dfSub$peak_id[i] = paste0("peak_", i)
  }
  dfSub$V5 = 1
  dfSub$strand = "+"
  write.table(dfSub, outFile, col.names = F, row.names = F, quote = F, sep = "\t")
  
}


# need conda activate homer
findMotifs = function(inFile, outDir) {
  dir.create(outDir, recursive = T, showWarnings = F)
  cmd_str = paste0("findMotifsGenome.pl ", inFile, " mm10 ", outDir, "/ -size 500 -mask")
  system(cmd_str)
}

#makeBed(inFile = "../chd_chd7_shun.txt", outFile = "../chd_chd7_shun.bed", header =F)
findMotifs(inFile = "../chd_chd7_shun.bed" , outDir = "motifs/homer_motif_chd7PeaksShun")

#makeBed(inFile = "../h3k_chd7_shun.txt", outFile = "../h3k_chd7_shun.bed", header =F)
#findMotifs(inFile = "../h3k_chd7_shun.bed" , outDir = "motifs/homer_motif_H3KPeaksShun")

#makeBed(inFile = "../chd_all_shun.txt", outFile = "../chd_all_shun.bed", header =T)
#makeBed(inFile = "../h3k_all_shun.txt", outFile = "../h3k_all_shun.bed", header =T)