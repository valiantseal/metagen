library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(BSgenome.Mmusculus.UCSC.mm10)
library(TFBSTools)
library(JASPAR2022)
library(karyoploteR)
library(monaLisa)

setwd("/home/flyhunter/Kai/Chipseq/pnas/process_custom")


startR = c(8647705, 8666955)
endR = c(8648095, 8667345)

curStart = startR[1]
curEnd = endR[1]

regStart = curStart
regEnd = curEnd

curChr = "chr4"

curRange = paste0(curChr, ":", regStart, "-", regEnd)
gr <- toGRanges(curRange)
seqs <- getSeq(BSgenome.Mmusculus.UCSC.mm10, gr)



####

#pfms <- getMatrixByID(JASPAR2022, c("MA0143.1", "MA0143.2", "MA0143.3"))
pfms <- getMatrixByID(JASPAR2022, c("MA0087.1", "MA0087.2"))

pwms <- toPWM(pfms)

res <- findMotifHits(query = pwms,
                     subject = seqs,
                     min.score = 7,
                     method = "matchPWM",
                     BPPARAM = BiocParallel::SerialParam())

#unique(res$pwmname)
res_df = data.frame(res)
length(unique(res_df$pwmname))
length(unique(res_df$start))

input_string <- as.character(seqs)

input_string = tolower(input_string)


showRange = function(imput_string, motifRes) {
  startPos = unique(res_df$start)
  endPos = unique(res_df$end)
  curChars = strsplit(input_string, "")[[1]]
  for (j in 1:length(startPos)) {
    curStartPos = startPos[j]
    curEndPos = endPos[j]
    for ( i in curStartPos:curEndPos) {
      curChars[i] = toupper(curChars[i])
    }
    resStr = paste(curChars, collapse = "")
  }

  return(resStr)
}

sox2Range = showRange(imput_string=imput_string, motifRes=res_df)

sox5Range = showRange(imput_string=imput_string, motifRes=res_df)

write.table(sox2Range, paste0(curChr, "_", regStart, "-", regEnd, "_Sox2Motif.txt"), col.names = F, row.names = F, quote = F)
write.table(sox5Range, paste0(curChr, "_", regStart, "-", regEnd, "_Sox5Motif.txt"), col.names = F, row.names = F, quote = F)
