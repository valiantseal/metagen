library(foreach)
library(doParallel)
library(readr)
library(dplyr)
library(stringr)

useCores <- 15

blastNames <- c('staxids', 'qseqid', 'sseqid', 'stitle', 'pident', 'length', 'mismatch',
                'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore')

filterReads <- function(df) {
  x <- readr::read_delim(df, delim = "\t", escape_double = FALSE, col_names = FALSE, comment = "#", trim_ws = TRUE)
  if (nrow(x) > 0) {
    colnames(x) <- blastNames
    x <- x[!is.na(x$bitscore), ]

    x <- x %>%
      filter(!grepl('Synthetic|clone', stitle, ignore.case = TRUE)) %>%
      arrange(qseqid, desc(bitscore)) %>%
      group_by(qseqid) %>%
      mutate(staxids = as.character(staxids)) %>%
      summarize(BID1 = first(staxids),
                BID2 = nth(unique(staxids), 2, default = NA_character_),
                BID3 = nth(unique(staxids), 3, default = NA_character_)) %>%
      ungroup()

    x <- x %>%
      mutate(BID1 = if_else(str_detect(BID1, ";"),
                            str_extract(BID1, "^[^;]+"),
                            BID1),
             BID2 = if_else(BID2 == BID1, NA_character_, BID2),
             BID3 = if_else(BID3 == BID1 | BID3 == BID2, NA_character_, BID3))
    return(x)
  } else {
    return(NULL)
  }
}

allDir <- list.files('./process')

for (sampleDir in allDir) {
  sampleDirPath <- paste0('./process/', sampleDir, '/')
  setwd(sampleDirPath)
  cl <- makeCluster(useCores, type = "FORK")
  registerDoParallel(cl)

  targDir <- './splitSeq10K/'
  filesList <- list.files(targDir)
  targSamples <- paste0(targDir, filesList, '/NtV4_blast.results')
  sampleName <- read.table('sample.name', F)

  blastComb <- foreach(i = targSamples, .combine = rbind, .packages = c('readr', 'dplyr', 'stringr')) %dopar% {
    filterReads(i)
  }

  if (!is.null(blastComb)) {
    blastComb <- blastComb %>%
      rename(Read = qseqid) %>%
      select(Read, BID1, BID2, BID3) %>%
      mutate(Sample = sampleName$V1)

    outFile <- 'blastResTop_3.tsv'
    write.table(blastComb, file = outFile, col.names = T, row.names = F, quote = T, sep = '\t')
  }

  stopCluster(cl)
  print(paste0(sampleDirPath, '________done!'))
  setwd('../../')
}
