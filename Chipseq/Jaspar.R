library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(BSgenome.Mmusculus.UCSC.mm10)
library(TFBSTools)
library(JASPAR2022)
library(karyoploteR)
library(monaLisa)

opts <- list()

opts[["species"]] <- 10090
#opts[["tax_group"]] <- "vertebrates"

PFMatrixList <- getMatrixSet(JASPAR2022, opts)

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

pwms <- toPWM(PFMatrixList)
#pwms


#name(pwms)

res <- findMotifHits(query = pwms,
                     subject = seqs,
                     min.score = 7.0,
                     method = "matchPWM",
                     BPPARAM = BiocParallel::SerialParam())

#unique(res$pwmname)

res_df = data.frame(res)
length(unique(res_df$pwmname))

write.csv(res_df, paste0("motifs/motifs_chr4_", curStart, "_", curEnd, "_jaspar2022_mouse.csv"), row.names = F)