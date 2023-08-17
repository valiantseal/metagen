library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(BSgenome.Mmusculus.UCSC.mm10)
library(TFBSTools)
library(JASPAR2022)
library(karyoploteR)
library(monaLisa)

setwd("/home/flyhunter/Kai/Chipseq/pnas/process_custom")

opts <- list()

#opts[["species"]] <- 10090
opts[["tax_group"]] <- "vertebrates"

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
# originally used min score 7
res <- findMotifHits(query = pwms,
                     subject = seqs,
                     min.score = "80%",
                     method = "matchPWM",
                     BPPARAM = BiocParallel::SerialParam())

#unique(res$pwmname)

res_df = data.frame(res)
length(unique(res_df$pwmname))

write.csv(res_df, paste0("motifs/motifs_chr4_", curStart, "_", curEnd, "_jaspar2022_80%_vertebrates.csv"), row.names = F)


## custom motifs
pfms <- getMatrixByID(JASPAR2022, c("MA1621.1", "MA1116.1"))
pwms <- toPWM(pfms)

res <- findMotifHits(query = pwms,
                     subject = seqs,
                     min.score = "80%",
                     method = "matchPWM",
                     BPPARAM = BiocParallel::SerialParam())

#unique(res$pwmname)
res_df = data.frame(res)
length(unique(res_df$pwmname))

# try with regular biostrins
startR = c(8647705, 8666955)
endR = c(8648095, 8667345)
curStart = startR[2]
curEnd = endR[2]
regStart = curStart
regEnd = curEnd
curChr = "chr4"
curRange = paste0(curChr, ":", regStart, "-", regEnd)
gr <- toGRanges(curRange)
motifs = c("MA1621.1", "MA1116.1")
pfms <- getMatrixByID(JASPAR2022, motifs[1])
pwms <- toPWM(pfms)
seqs <- getSeq(BSgenome.Mmusculus.UCSC.mm10, gr)
dna_seq = paste(seqs, collapse = "")
pwmMat = pwms@profileMatrix
matches <- Biostrings::matchPWM(pwmMat, dna_seq, min.score="10%")
matches

### try with homer
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

pfms <- getMatrixByID(JASPAR2022, c("MA1621.1", "MA1116.1"))
pwms <- toPWM(pfms)


res <- findMotifHits(query = pwms,
                     subject = seqs,
                     method = "homer2",
                     homerfile = "/home/flyhunter/miniconda3/envs/homer/bin/homer",
                     min.score = 6)

#unique(res$pwmname)

res_df = data.frame(res)
length(unique(res_df$pwmname))
